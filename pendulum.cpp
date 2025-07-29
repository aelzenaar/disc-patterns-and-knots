#define _GLIBCXX_DEBUG
#include <complex>
using namespace std::complex_literals;
#include <numbers>

#include <iostream>
#include <array>
#include <iomanip>
#include <thread>

#include <eigen3/Eigen/Dense>

using RR = double;
using CC = std::complex<RR>;
using Mob = Eigen::Matrix<CC, 2,2>;



inline std::array<Mob, 6> representation(RR x, RR y)
{
  const CC Eix = std::exp(1.i * x * std::numbers::pi_v<RR>);
  const CC Eiy = std::exp(1.i * y * std::numbers::pi_v<RR>);
  auto m1 = Mob{{-Eix, 2.*Eix*Eix}, {-2., 3.*Eix}};
  auto m2 = Mob{{-4.*Eix - Eiy, 8.*Eix*Eix + 8.*Eix*Eiy + 2.*Eiy*Eiy}, {- 2., 4.*Eix + 3.*Eiy}};
  auto m3 = Mob{{- 7.*Eix*Eiy - 4.*Eiy*Eiy - 4.*Eix*Eix, 2.*Eix*Eix*Eiy + 2.*Eix*Eiy*Eiy}, {- 2.*Eix - 2.*Eiy, Eix*Eiy}};
  m1 = m1/std::sqrt(m1.determinant());
  m2 = m2/std::sqrt(m2.determinant());
  m3 = m3/std::sqrt(m3.determinant());

  return {
    m1, m2, m3,
    m1.inverse(), m2.inverse(), m3.inverse()
  };
}

using IndexList = std::vector<size_t>;
std::array<IndexList, 6> allowed_followers = {
  IndexList{0,1,2,4,5},
  IndexList{0,1,2,3,5},
  IndexList{0,1,2,3,4},
  IndexList{1,2,3,4,5},
  IndexList{0,2,3,4,5},
  IndexList{0,1,3,4,5},
};


const size_t depth = 9;

using WordTree = std::vector<std::vector<std::pair<size_t, Mob> > >;

WordTree word_tree_from_matrices(std::array<Mob, 6> matrices)
{
  WordTree base_point_list;
  {
    std::vector<std::pair<size_t, Mob> > list0;
    for(auto i = 0; i < matrices.size(); i++)
        list0.push_back(std::make_pair(i, matrices[i]));
    base_point_list.push_back(list0);
  }

  for(auto height = 1 ; height <= depth; height++)
  {
    std::vector<std::pair<size_t, Mob> > list0;
    for(auto& previous_word : base_point_list.back())
    {
      for(auto n : allowed_followers[previous_word.first])
        list0.push_back(std::make_pair(n, matrices[n] * previous_word.second));
    }
    base_point_list.push_back(list0);
  }

  return base_point_list;
}

inline RR real_length(Mob f)
{
  auto inp = f.trace();
  return std::abs(2*log(std::abs((inp + std::sqrt(inp*inp-4.))/2.)));
}

std::vector<RR> lengths_from_word_tree(WordTree matrices)
{
  std::vector<RR> lengths;
  for(auto& x : matrices)
    for(auto& y : x)
      lengths.push_back(real_length(y.second));
  return lengths;
}

void vertical_strip(const std::vector<RR>& base_point_lengths,
                    unsigned short px_width,
                    unsigned short px_height,
                    unsigned short start_px_x,
                    unsigned short end_px_x,
                    CC view_corner_bl,
                    CC view_corner_tr)
{
  for(unsigned short px_x = start_px_x; px_x < end_px_x; px_x++)
  {
    for(unsigned short px_y = 0; px_y < px_height; px_y++)
    {
      std::cerr << "Computing heights, x="<<px_x<<" y="<<px_y<<", approx perc = "<< std::fixed << std::setprecision(2) << (RR)(px_x*px_height + px_y)/(RR)(px_height*px_width) << std::endl;
      RR x = ((RR)px_x/(RR)px_width * std::abs(view_corner_tr.real() - view_corner_bl.real()) + view_corner_bl.real());
      RR y = ((RR)px_y/(RR)px_height * std::abs(view_corner_tr.imag() - view_corner_bl.imag()) + view_corner_bl.imag());

      auto pixel_lengths = lengths_from_word_tree(word_tree_from_matrices(representation(x,y)));

      RR shortest = std::numeric_limits<RR>::infinity();
      for(size_t i = 0; i < base_point_lengths.size(); i++)
      {
        if(base_point_lengths[i] < 1e-5) continue;
        if(pixel_lengths[i] < shortest) shortest = pixel_lengths[i];
      }
      std::printf("%f,%f,%f\n",x,y,shortest);
      // std::stringstream stream;
      // stream << x << "," << y << "," << std::fixed << std::setprecision(6) << shortest << std::endl;
      // std::cout<<stream.str();
    }
  }
}


int main()
{
  auto tree = word_tree_from_matrices(representation(0,0));
  auto base_point_lengths = lengths_from_word_tree(tree);

  for(auto list : tree)
    std::cerr<<list.size()<<" ";
  std::cerr<<std::endl;


  unsigned short px_width = 400;
  unsigned short px_height = 400;

  std::vector<std::thread> threads;
  for(unsigned short lr = 0; lr < 40; lr++)
    threads.push_back(std::thread([=, &base_point_lengths]() {vertical_strip(base_point_lengths, px_width, px_height, lr*10, lr*10 + 10, 0., 2.+2.i);} ));
  for(auto& t:threads) t.join();

  return 0;
}
