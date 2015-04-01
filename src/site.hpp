/* @(#)site.hpp
 */

#ifndef _SITE_HPP
#define _SITE_HPP 1

#include <vector>
#include <string>
#include <cassert>

class Site_base {
public:
  std::string chr;
  uint32_t pos = 0;
  std::vector<std::string> all;
  std::string id;

  void init(std::string chr, uint32_t pos, std::vector<std::string> alls,
            std::string id = ".") {
    assert(alls.size() == 2);
    this->chr = std::move(chr);
    this->pos = pos;
    this->id = std::move(id);
    all = std::move(alls);
  }

  bool operator==(const Site_base &lhs) {
    if (chr != lhs.chr || pos != lhs.pos)
      return false;
    if (all.size() != lhs.all.size())
      return false;
    for (size_t i = 0; i < all.size(); ++i)
      if (all[i] != lhs.all[i])
        return false;
    if (id != lhs.id)
      return false;
    return true;
  }
  bool operator!=(const Site_base &lhs) { return !(this->operator==(lhs)); }

  bool empty() {
    if (chr.empty() || all.empty())
      return true;
    return false;
  }

  void flipStrand() { std::swap(all[0], all[1]); }
};

class Site : public Site_base {
public:
  std::vector<double> hap;
  double weight;

  void flipStrand() {
    Site_base::flipStrand();
    for (auto &h : hap)
      h = -(h - 0.5) + 0.5;
  }
};

#endif /* _SITE_HPP */
