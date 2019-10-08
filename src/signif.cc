#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <array>

#include <nlohmann/json.hpp>

#include "ivanp/string.hh"
#include "ivanp/wls.hh"
#include "ivanp/math/polynomial.hh"

#define TEST(var) \
  cout << "\033[36m" #var "\033[0m = " << (var) << endl;

using namespace std;
using ivanp::cat;

size_t utf8len(const char* s) {
  size_t len = 0;
  while (*s) len += (*s++ & 0xc0) != 0x80;
  return len;
}
size_t utf8len(const std::string& s) { return utf8len(s.c_str()); }

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cout << "usage: " << argv[0] << " in.json\n";
    return 1;
  }
  nlohmann::json in;
  ifstream(argv[1]) >> in;

  const double lumi = in["lumi"];
  cout << "Lumi = " << lumi << endl;

  vector<tuple<string,vector<double>>> vars = in["bins"];
  vector<unsigned> nbins; nbins.reserve(vars.size());
  for (const auto& var : vars) {
    const auto& edges = get<1>(var);
    nbins.push_back(edges.size()-1);
    cout << get<0>(var) << ": " << nbins.back() << ' '
      << edges.front() << " - " << edges.back();
    if (edges.size()<20) {
      cout << '\n';
      for (auto edge : edges) cout << ' ' << edge;
    }
    cout << endl;
  }

  vector<double> xs; xs.reserve(nbins[0]);
  const auto& edges = get<1>(vars[0]);
  auto skip = [&](unsigned i){ // skip signal region
    return !(edges[i]<121 || 129<edges[i+1]);
  };
  for (unsigned i=0; i<nbins[0]; ++i) {
    if (skip(i)) continue;
    xs.push_back((edges[i]+edges[i+1])/2);
  }
  // TEST(xs.size())
  cout << endl;

  vector<double> cs(3), tcs(cs.size());
  vector<double> A;
  A.reserve(cs.size()*xs.size());
  for (unsigned p=0; p<cs.size(); ++p)
    for (auto x : xs)
      A.push_back(pow(x,p));

  vector<vector<string>> results;
  results.emplace_back();
  for (unsigned i=vars.size()-1; i; --i)
    results[0].emplace_back(get<0>(vars[i]));
  for (const auto& str : {"sig","bkg","s/√(s+b)"})
    results[0].emplace_back(str);

  vector<array<double,2>> data = in["data"], mc = in["mc"];
  vector<double> ys(xs.size()), us(xs.size());
  vector<unsigned> ii(nbins.size());
  double sig = 0;
  for (unsigned i=0, j=0, k=0; i<data.size(); ++i) {
    if (skip(j++)) {
      sig += mc[i][0];
      continue;
    }
    const double y = data[i][0];
    ys[k] = log(y);
    us[k] = data[i][1]/y;
    if (j == nbins[0]) {
      results.emplace_back();
      auto res = [&r=results.back()](const auto&... args){
        r.emplace_back(cat(args...));
      };

      for (unsigned v=nbins.size()-1; v; --v) {
        const auto& edges = get<1>(vars[v]);
        res("[",edges[ii[v]],',',edges[ii[v]+1],")");
      }
      ++ii[1];
      for (unsigned v=1; v<nbins.size(); ++v) {
        if (ii[v]==nbins[v]) {
          ii[v] = 0;
          ++ii[v+1];
        }
      }

      sig *= lumi;
      // TEST(sig)
      res(sig);

      ivanp::wls(A.data(),ys.data(),us.data(),xs.size(),cs.size(),cs.data());
      // for (auto c : cs)
      //   TEST(c)

      ivanp::math::poly::transform_coords(121,8,3,cs.data(),tcs.data());

      const double
        bkg = (std::exp(tcs[0])/3)
            * ( tcs[1]*(tcs[1]*(tcs[1]+4) + 6*(tcs[2]+2)) + 8*(tcs[2]+3) );
      // TEST(bkg)
      res(bkg);

      const double signif = sig/sqrt(sig+bkg);
      // TEST(signif)
      res(signif);

      sig = 0;
      j = 0;
      k = 0;
    } else ++k;
  }

  vector<unsigned> width(results[0].size());
  for (const auto& r : results) {
    for (unsigned i=0; i<r.size(); ++i) {
      const auto len = utf8len(r[i]);
      if (width[i] < len) width[i] = len;
    }
  }
  for (unsigned j=0; j<results.size(); ++j) {
    for (unsigned i=0; i<width.size(); ++i) {
      if (i) cout << "  ";
      if (i<vars.size()-1) cout << left;
      else cout << right;
      cout << setw(width[i]) << results[j][i];
    }
    cout << '\n';
    if (!(j%nbins[1])) cout << '\n';
  }
  cout << endl;
}
