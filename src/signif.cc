#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <array>

#include <nlohmann/json.hpp>

#include "ivanp/wls.hh"
#include "ivanp/math/polynomial.hh"

#define TEST(var) \
  cout << "\033[36m" #var "\033[0m = " << (var) << endl;

using namespace std;

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
  TEST(xs.size())

  vector<double> cs(3), tcs(cs.size());
  vector<double> A;
  A.reserve(cs.size()*xs.size());
  for (unsigned p=0; p<cs.size(); ++p)
    for (auto x : xs)
      A.push_back(pow(x,p));

  vector<array<double,2>> data = in["data"], mc = in["mc"];
  vector<double> ys(xs.size()), us(xs.size());
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
      ivanp::wls(A.data(),ys.data(),us.data(),xs.size(),cs.size(),cs.data());
      for (auto c : cs)
        TEST(c)

      ivanp::math::poly::transform_coords(121,8,3,cs.data(),tcs.data());

      const double
        bkg = (std::exp(tcs[0])/3)
            * ( tcs[1]*(tcs[1]*(tcs[1]+4) + 6*(tcs[2]+2)) + 8*(tcs[2]+3) );
      TEST(bkg)

      sig *= lumi;
      TEST(sig)

      const double signif = sig/sqrt(sig+bkg);
      TEST(signif)

      sig = 0;
      j = 0;
      k = 0;
    } else ++k;
  }
}
