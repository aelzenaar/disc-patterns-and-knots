// #define _GLIBCXX_DEBUG
#include <complex>
using namespace std::complex_literals;

#include <iostream>
#include <array>
#include <iomanip>
#include <random>
#include <thread>

#include "beta_dist.hpp"

#include <eigen3/Eigen/Dense>

using RR = long double;
using CC = std::complex<RR>;
using Mob = Eigen::Matrix<CC, 2,2>;

inline std::array<Mob, 4> representation(CC tX, CC tY, CC tXY)
{
  const CC v = std::sqrt(4.L-tXY*tXY);
  return {
    Mob{{.5L*tX + .5il*tXY, - .5L*tX - .5L*v}, {-.5L*tX + .5L*v, .5L*tX - .5il*tXY}},
    Mob{{.5L*tY - 1.il, .5L*tY}, {.5L*tY, .5L*tY + 1.il}},
    Mob{{.5L*tX - .5il*tXY, .5L*tX + .5L*v}, {.5L*tX - .5L*v, .5L*tX + .5il*tXY}},
    Mob{{.5L*tY + 1.il, -.5L*tY}, {-.5L*tY, .5L*tY - 1.il}}
  };
}


std::discrete_distribution<> any_follower({0,1,2,3});

using IndexList = std::vector<size_t>;
std::array<std::discrete_distribution<>, 4> allowed_followers = {
  std::discrete_distribution<>({0,1,3}), // X can be followed by XYy
  std::discrete_distribution<>({0,1,2}), // Y can be followed by XYx
  std::discrete_distribution<>({1,2,3}), // x can be followed by Yxy
  std::discrete_distribution<>({0,2,3}) // y can be followed by Xxy
};


const size_t min_length = 1;
const size_t max_length = 20;
const size_t count = 20'000'000;

using WordTree = std::vector<std::pair<size_t, Mob> >;


std::vector<IndexList> random_words()
{
  std::vector<IndexList> list;

  std::random_device rd;
  std::mt19937 gen(rd());
  beta_distribution<RR> beta(5, 4);
  for(auto i = 0; i < count; ++i)
  {
    auto len = min_length + beta(gen) * (max_length - min_length);
    IndexList word;
    for(auto j = 0; j < len; ++j)
    {
      if(word.size() == 0) word.push_back(any_follower(gen));
      else word.push_back(allowed_followers[word.back()](gen));
    }
    list.push_back(word);
  }
  return list;
}


WordTree word_tree_from_matrices(std::array<Mob, 4> matrices, std::vector<IndexList> word_list)
{
  WordTree base_point_list;
  for(auto word : word_list)
  {
    Mob m = Mob::Identity();
    for(auto w : word) m *= matrices[w];
    base_point_list.push_back(std::make_pair(word.size(), m));
  }

  return base_point_list;
}

inline RR real_length(Mob f)
{
  auto inp = f.trace();
  return std::abs(2*log(std::abs((inp + std::sqrt(inp*inp-4.L))/2.L)));
}

std::vector<RR> lengths_from_word_tree(WordTree matrices)
{
  std::vector<RR> lengths;
  for(auto& y : matrices)
      lengths.push_back(real_length(y.second));
  return lengths;
}

void vertical_strip(const std::vector<RR>& base_point_lengths,
                    std::vector<IndexList> word_list,
                    unsigned short px_width,
                    unsigned short px_height,
                    unsigned short start_px_x,
                    unsigned short end_px_x,
                    CC view_corner_bl,
                    CC view_corner_tr)
{
  std::stringstream filename;
  filename << "schottky_slice_hard_rand_out/x" << start_px_x << "_" << end_px_x << ".csv" ;
  std::FILE* output_file = std::fopen(filename.str().c_str(), "w");
  for(unsigned short px_x = start_px_x; px_x < end_px_x; px_x++)
  {
    for(unsigned short px_y = 0; px_y < px_height; px_y++)
    {
      std::cerr << "Computing heights, x="<<px_x<<" y="<<px_y<<", approx perc = "<< std::fixed << std::setprecision(2) << (RR)(px_x*px_height + px_y)/(RR)(px_height*px_width) << std::endl;
      RR x = ((RR)px_x/(RR)px_width * std::abs(view_corner_tr.real() - view_corner_bl.real()) + view_corner_bl.real());
      RR y = ((RR)px_y/(RR)px_height * std::abs(view_corner_tr.imag() - view_corner_bl.imag()) + view_corner_bl.imag());

      const CC t = x+1il*y;

      auto pixel_lengths = lengths_from_word_tree(word_tree_from_matrices(representation(t*(0.76069L + 0.857874il)+(1.L-t)*(2.L- 1.il),
                                                                                         t*(-0.76069L - 0.857874il)+(1.L-t)*(- 1.il),
                                                                                         t*(-3.05898L - 1.9751il)+(1.L-t)*(-2.L-2.il)),
                                                                          word_list));

      RR shortest = std::numeric_limits<RR>::infinity();
      for(size_t i = 0; i < base_point_lengths.size(); i++)
      {
        if(base_point_lengths[i] < 1e-5) continue;
        if(pixel_lengths[i] < shortest) shortest = pixel_lengths[i];
      }
      std::fprintf(output_file,"%Lf,%Lf,%Lf\n",x,y,shortest);
    }
  }
  std::fclose(output_file);
}


int main()
{
  std::vector<IndexList> word_list = random_words();
  auto tree = word_tree_from_matrices(representation(-3,-3,-3), word_list);
  auto base_point_lengths = lengths_from_word_tree(tree);

  unsigned short n_threads = 80;
  unsigned short lines_per_thread = 16;
  CC bottom_left = -2.L-2.il;
  CC top_right = 2.L+2.il;

  unsigned short px_width = n_threads*lines_per_thread;
  unsigned short px_height = px_width;

  std::vector<std::thread> threads;
  for(unsigned short lr = 0; lr < n_threads; lr++)
    threads.push_back(std::thread([=, &base_point_lengths]() {vertical_strip(base_point_lengths, word_list, px_width, px_height, lr*lines_per_thread, lr*lines_per_thread + lines_per_thread, bottom_left, top_right);} ));
  for(auto& t:threads) t.join();

  return 0;
}
