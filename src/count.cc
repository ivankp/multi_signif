#include <iostream>

#include "ivanp/io/mem_file.hh"
#include "ivanp/timed_counter.hh"
#include "reader2.hh"

#define TEST(var) \
  std::cout << "\033[36m" #var "\033[0m = " << var << std::endl;

using std::cout;
using std::endl;
using std::cerr;

using vec4 = ivanp::vec4<double>;

float lumi=0, weight=1;
uint32_t nevents_total = 0;
uint8_t njets = 0;

int main(int argc, char* argv[]) {
  reader read(argv[1]);

  char is_mc;
  read(is_mc);
  if (is_mc!='m' && is_mc!='d') {
    cerr << "file begins with \'" << is_mc << "\'" << endl;
    return 1;
  }
  is_mc = (is_mc=='m');

  if (!is_mc) {
    read(lumi);
    TEST(lumi);
    weight = 1;
  }
  read(nevents_total);
  TEST(nevents_total);

  size_t left_side=0, right_side=0;

  vec4 y[2], yy;
  { ivanp::timed_counter<> ent;
    for (; read; ++ent) {
      if (is_mc) read(weight);
      read(y[0]);
      read(y[1]);
      yy = y[0] + y[1];
      auto myy = yy.m();

      if (myy < 121) ++left_side;
      else if (myy > 129) ++right_side;

      read(njets);
      if (njets>4) njets = 4;
      read.skip(njets*sizeof(float[4]));
    }
  }

  cout << endl;
  TEST(left_side)
  TEST(right_side)
  TEST(left_side+right_side)
}
