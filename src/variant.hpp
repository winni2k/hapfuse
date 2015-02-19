/* @(#)variant.hpp
 */

#ifndef _VARIANT_HPP
#define _VARIANT_HPP 1

class Variant {

private:
  std::string m_chrom;
  std::string m_id;
  std::string m_ref;
  std::string m_alt;
  unsigned m_pos;

  Variant(std::string chrom, std::string id, unsigned pos, std::string ref,
          std::string alt)
      : m_chrom(std::move(chrom)), m_id(std::move(id)), m_pos(pos),
        m_ref(std::move(ref)), m_alt(std::move(alt)) {}

    
}

#endif /* _VARIANT_HPP */
