#define _GLIBCXX_DEBUG
#include <complex>
using namespace std::complex_literals;

#include <iostream>
#include <array>
#include <iomanip>
#include <thread>

#include <eigen3/Eigen/Dense>

using RR = long double;
using CC = std::complex<RR>;
using Mob = Eigen::Matrix<CC, 2,2>;

inline std::array<Mob, 6> representation(CC alpha, CC beta, CC lambda)
{
  return {
    Mob{{1.0L, alpha}, {0, 1.0L}},
    Mob{{1.0L, beta}, {0, 1.0L}},
    Mob{{lambda, lambda*lambda-1.0L}, {1.0L, lambda}},
    Mob{{1.0L, -alpha}, {0, 1.0L}},
    Mob{{1.0L, -beta}, {0, 1.0L}},
    Mob{{lambda, -lambda*lambda+1.0L}, {-1.0L, lambda}}
  };
}

using IndexList = std::vector<size_t>;
std::array<IndexList, 6> allowed_followers = {
  IndexList{0,1,2,4,5}, // P can be followed by PQM qm
  IndexList{1,2,5}, // Q can be followed by QMm
  IndexList{0,1,2,3,4}, // M can be followed by PQMpq
  IndexList{1,2,3,4,5}, // p can be followed by QMpqm
  IndexList{2,4,5}, // q can be followed by Mqm
  IndexList{0,1,3,4,5} // m can be followed by PQpqm
};


const size_t depth = 10;

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
  return std::abs(2*log(std::abs((inp + std::sqrt(inp*inp-4.L))/2.L)));
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
                    const unsigned short px_width,
                    const unsigned short px_height,
                    const unsigned short start_px_x,
                    const unsigned short end_px_x,
                    const CC view_corner_bl,
                    const CC view_corner_tr)
{
  for(unsigned short px_x = start_px_x; px_x < end_px_x; px_x++)
  {
    for(unsigned short px_y = 0; px_y < px_height; px_y++)
    {
      std::cerr << "Computing heights, x="<<px_x<<" y="<<px_y<<", approx perc = "<< std::fixed << std::setprecision(2) << (RR)(px_x*px_height + px_y)/(RR)(px_height*px_width) << std::endl;
      const RR x = ((RR)px_x/(RR)px_width * std::abs(view_corner_tr.real() - view_corner_bl.real()) + view_corner_bl.real());
      const RR y = ((RR)px_y/(RR)px_height * std::abs(view_corner_tr.imag() - view_corner_bl.imag()) + view_corner_bl.imag());

      const CC t = x+1il*y;

//       In[15]:= Solve[{Tr[M]==2, M.Inverse[P].M.Inverse[Q] == {{-1,0},{0,-1}}, Tr[Inverse[P].M.Inverse[Q]] == -2}, {a,b,l}]
//
//       Out[15]= {{a -> 2, b -> 2, l -> 1}}
//


      auto pixel_lengths = lengths_from_word_tree(word_tree_from_matrices(representation(2.l+ 2.il*t,
                                                                                         2.l-2.il*t,
                                                                                         1.l)));

      RR shortest = std::numeric_limits<RR>::infinity();
      for(size_t i = 0; i < base_point_lengths.size(); i++)
      {
        if(base_point_lengths[i] < 1e-5) continue;
        if(pixel_lengths[i] < shortest) shortest = pixel_lengths[i];
      }
      std::printf("%Lf,%Lf,%Lf\n",x,y,shortest);
      // std::stringstream stream;
      // stream << x << "," << y << "," << std::fixed << std::setprecision(6) << shortest << std::endl;
      // std::cout<<stream.str();
    }
  }
}


int main()
{
  auto tree = word_tree_from_matrices(representation(2.l + 2.il*1.5l, 2.l - 2.il*1.5l, 1.l));
  auto base_point_lengths = lengths_from_word_tree(tree);

  for(auto list : tree)
    std::cerr<<list.size()<<" ";
  std::cerr<<std::endl;


  unsigned short px_width = 400;
  unsigned short px_height = 400;

  std::vector<std::thread> threads;
  for(unsigned short lr = 0; lr < 40; lr++)
    threads.push_back(std::thread([=, &base_point_lengths]() {vertical_strip(base_point_lengths, px_width, px_height, lr*10, lr*10 + 10, -1.L-2.il, 3.L+2.il);} ));
  for(auto& t:threads) t.join();

  return 0;
}
