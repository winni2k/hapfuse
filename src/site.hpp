/* @(#)site.hpp
 */

#ifndef _SITE_HPP
#define _SITE_HPP 1

#include <vector>
#include <string>


class Site_base {
public:
  std::string chr;
  uint32_t pos = 0;
  std::vector<std::string> all;

  void init(std::string chr, uint32_t pos, std::vector<std::string> alls) {
    this->chr = std::move(chr);
    this->pos = pos;
    all = std::move(alls);
    //    for (auto &a : alls)
    //      all.push_back(std::move(a));
  }

  bool operator==(const Site_base &lhs) {
    if (chr != lhs.chr || pos != lhs.pos)
      return false;
    if (all.size() != lhs.all.size())
      return false;
    for (size_t i = 0; i < all.size(); ++i)
      if (all[i] != lhs.all[i])
        return false;
    return true;
  }
  bool operator!=(const Site_base &lhs) { return !(this->operator==(lhs)); }

  bool empty() {
    if (chr.empty() || all.empty())
      return true;
    return false;
  }

};

class Site: public Site_base {
public:
  std::vector<double> hap;
  double weight;

  void flipStrand(){
    std::swap(all[0],all[1]);
    for(auto &h : hap)
      h = -(h - 0.5) + 0.5;
  }
};

#endif /* _SITE_HPP */
