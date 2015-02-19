/* @(#)site.hpp
 */

#ifndef _SITE_HPP
#define _SITE_HPP 1

class Site {
public:
  std::vector<double> hap;
  std::string chr;
  uint32_t pos = 0;
  double weight;
  std::vector<std::string> all;

  void init(std::string chr, uint32_t pos,
            std::vector<std::string> alls) {
    this->chr = std::move(chr);
    this->pos = pos;
    for (auto &a : alls)
      all.push_back(std::move(a));
  }

  bool operator==(const Site &lhs) {
    if (chr != lhs.chr || pos != lhs.pos)
      return false;
    if (all.size() != lhs.all.size())
      return false;
    for (size_t i = 0; i < all.size(); ++i)
      if (all[i] != lhs.all[i])
        return false;
    return true;
  }

  bool empty() {
    if (chr.empty() || all.empty())
      return true;
    return false;
  }
};

#endif /* _SITE_HPP */
