/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2020 Lucas Czech and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief Implementation of bitvector functions.
 *
 * @file
 * @ingroup utils
 */

#include "bitvector.hpp"

#include <stdexcept>

namespace genesis {
namespace utils {

// =============================================================================
//     Constants
// =============================================================================

const Bitvector::IntType Bitvector::all_0_ = 0ul;
const Bitvector::IntType Bitvector::all_1_ = (((1ul << 32) - 1) << 32)  +  ((1ul << 32) - 1);

const std::array<Bitvector::IntType, Bitvector::IntSize> Bitvector::bit_mask_ =
{{
    1ul << 0,   1ul << 1,   1ul << 2,   1ul << 3,   1ul << 4,   1ul << 5,   1ul << 6,   1ul << 7,
    1ul << 8,   1ul << 9,   1ul << 10,  1ul << 11,  1ul << 12,  1ul << 13,  1ul << 14,  1ul << 15,
    1ul << 16,  1ul << 17,  1ul << 18,  1ul << 19,  1ul << 20,  1ul << 21,  1ul << 22,  1ul << 23,
    1ul << 24,  1ul << 25,  1ul << 26,  1ul << 27,  1ul << 28,  1ul << 29,  1ul << 30,  1ul << 31,
    1ul << 32,  1ul << 33,  1ul << 34,  1ul << 35,  1ul << 36,  1ul << 37,  1ul << 38,  1ul << 39,
    1ul << 40,  1ul << 41,  1ul << 42,  1ul << 43,  1ul << 44,  1ul << 45,  1ul << 46,  1ul << 47,
    1ul << 48,  1ul << 49,  1ul << 50,  1ul << 51,  1ul << 52,  1ul << 53,  1ul << 54,  1ul << 55,
    1ul << 56,  1ul << 57,  1ul << 58,  1ul << 59,  1ul << 60,  1ul << 61,  1ul << 62,  1ul << 63
}};

const std::array<Bitvector::IntType, Bitvector::IntSize> Bitvector::ones_mask_ =
{{
    Bitvector::all_1_,       Bitvector::all_1_ >> 63, Bitvector::all_1_ >> 62, Bitvector::all_1_ >> 61,
    Bitvector::all_1_ >> 60, Bitvector::all_1_ >> 59, Bitvector::all_1_ >> 58, Bitvector::all_1_ >> 57,
    Bitvector::all_1_ >> 56, Bitvector::all_1_ >> 55, Bitvector::all_1_ >> 54, Bitvector::all_1_ >> 53,
    Bitvector::all_1_ >> 52, Bitvector::all_1_ >> 51, Bitvector::all_1_ >> 50, Bitvector::all_1_ >> 49,
    Bitvector::all_1_ >> 48, Bitvector::all_1_ >> 47, Bitvector::all_1_ >> 46, Bitvector::all_1_ >> 45,
    Bitvector::all_1_ >> 44, Bitvector::all_1_ >> 43, Bitvector::all_1_ >> 42, Bitvector::all_1_ >> 41,
    Bitvector::all_1_ >> 40, Bitvector::all_1_ >> 39, Bitvector::all_1_ >> 38, Bitvector::all_1_ >> 37,
    Bitvector::all_1_ >> 36, Bitvector::all_1_ >> 35, Bitvector::all_1_ >> 34, Bitvector::all_1_ >> 33,
    Bitvector::all_1_ >> 32, Bitvector::all_1_ >> 31, Bitvector::all_1_ >> 30, Bitvector::all_1_ >> 29,
    Bitvector::all_1_ >> 28, Bitvector::all_1_ >> 27, Bitvector::all_1_ >> 26, Bitvector::all_1_ >> 25,
    Bitvector::all_1_ >> 24, Bitvector::all_1_ >> 23, Bitvector::all_1_ >> 22, Bitvector::all_1_ >> 21,
    Bitvector::all_1_ >> 20, Bitvector::all_1_ >> 19, Bitvector::all_1_ >> 18, Bitvector::all_1_ >> 17,
    Bitvector::all_1_ >> 16, Bitvector::all_1_ >> 15, Bitvector::all_1_ >> 14, Bitvector::all_1_ >> 13,
    Bitvector::all_1_ >> 12, Bitvector::all_1_ >> 11, Bitvector::all_1_ >> 10, Bitvector::all_1_ >> 9,
    Bitvector::all_1_ >> 8,  Bitvector::all_1_ >> 7,  Bitvector::all_1_ >> 6,  Bitvector::all_1_ >> 5,
    Bitvector::all_1_ >> 4,  Bitvector::all_1_ >> 3,  Bitvector::all_1_ >> 2,  Bitvector::all_1_ >> 1
}};

const std::array<Bitvector::IntType, 4> Bitvector::count_mask_ =
{{
    0x5555555555555555,  //binary: 0101...
    0x3333333333333333,  //binary: 00110011...
    0x0f0f0f0f0f0f0f0f,  //binary: 4 zeros, 4 ones...
    0x0101010101010101   //the sum of 256 to the power of 0,1,2,3...
}};

// =============================================================================
//     Constructor and Rule of Five
// =============================================================================

Bitvector::Bitvector( std::string const& values )
    : Bitvector::Bitvector( values.size(), false )
{
    for( size_t i = 0; i < values.size(); ++i ) {
        switch( values[i] ) {
            case '0':
                break;
            case '1':
                set(i);
                break;
            default:
                throw std::invalid_argument(
                    "Cannot construct Bitvector from std::string that contains characters "
                    "other than 0 and 1."
                );
        }
    }
}

Bitvector::Bitvector( Bitvector const& other, size_t max_size )
{
    if( max_size > other.size() ) {
        max_size = other.size();
    }
    size_ = max_size;
    auto const ds = (size_ / IntSize) + (size_ % IntSize == 0 ? 0 : 1);
    assert( ds <= other.data_.size() );
    data_ = std::vector<IntType>( other.data_.begin(), other.data_.begin() + ds );
    unset_padding_();
}

// =============================================================================
//     Operators
// =============================================================================

Bitvector& Bitvector::operator &= (Bitvector const& rhs)
{
    if( size_ != rhs.size_ ) {
        throw std::runtime_error(
            "Cannot use operator `&` or `&=` on Bitvectors of different size. "
            "Use bitwise_and() instead."
        );
    }

    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] &= rhs.data_[i];
    }
    return *this;
}

bool Bitvector::operator < (const Bitvector &other) const
{
    if( size_ != other.size_ ) {
        throw std::runtime_error(
            "Cannot use operator `<` on Bitvectors of different size. "
        );
    }

    for (size_t i = 0; i < data_.size(); ++i) {
      if (data_[i] < other.data_[i]) {
        return true;
      } else if (data_[i] > other.data_[i]) {
        return false;
      }
    }
    return false;
}

Bitvector& Bitvector::operator |= (Bitvector const& rhs)
{
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] |= rhs.data_[i];
    }
    return *this;
}

Bitvector& Bitvector::operator ^= (Bitvector const& rhs)
{
    if( size_ != rhs.size_ ) {
        throw std::runtime_error(
            "Cannot use operator `^` or `^=` on Bitvectors of different size. "
            "Use bitwise_xor() instead."
        );
    }

    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] ^= rhs.data_[i];
    }
    return *this;
}

Bitvector Bitvector::operator ~ () const
{
    Bitvector cpy = Bitvector(*this);
    cpy.negate();
    return cpy;
}

Bitvector operator & (Bitvector const& lhs, Bitvector const& rhs)
{
    Bitvector result = Bitvector(lhs);
    result &= rhs;
    return result;
}

Bitvector operator | (Bitvector const& lhs, Bitvector const& rhs)
{
    Bitvector result = Bitvector(lhs);
    result |= rhs;
    return result;
}

Bitvector operator ^ (Bitvector const& lhs, Bitvector const& rhs)
{
    Bitvector result = Bitvector(lhs);
    result ^= rhs;
    return result;
}

bool Bitvector::operator == (const Bitvector &other) const
{
    if (size_ != other.size_) {
        return false;
    }
    for (size_t i = 0; i < data_.size(); ++i) {
        if (data_[i] != other.data_[i]) {
            return false;
        }
    }
    return true;
}

bool Bitvector::operator != (const Bitvector &other) const
{
    return !(*this == other);
}

// =============================================================================
//     Other Functions
// =============================================================================

size_t Bitvector::count() const
{
    size_t res = 0;
    for (IntType x : data_) {
        // put count of each 2 bits into those 2 bits
        x -= (x >> 1) & count_mask_[0];

        // put count of each 4 bits into those 4 bits
        x = (x & count_mask_[1]) + ((x >> 2) & count_mask_[1]);

        // put count of each 8 bits into those 8 bits
        x = (x + (x >> 4)) & count_mask_[2];

        // take left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
        res += (x * count_mask_[3]) >> 56;
    }

    // safe, but slow version...
    //~ size_t tmp = 0;
    //~ for (size_t i = 0; i < size_; ++i) {
        //~ if (get(i)) {
            //~ ++tmp;
        //~ }
    //~ }
    //~ assert(tmp == res);

    return res;
}


void Bitvector::negate()
{
    // flip all bits.
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] = ~ data_[i];
    }

    // reset the surplus bits at the end of the vector.
    unset_padding_();
}

void Bitvector::normalize()
{
    if (size_ > 0 && get(0)) {
        negate();
    }
}

void Bitvector::set_all( const bool value )
{
    // set according to flag.
    const auto v = value ? all_1_ : all_0_;
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] = v;
    }

    // if we initialized with true, we need to unset the surplus bits at the end!
    if (value) {
        unset_padding_();
    }
}

void Bitvector::unset_padding_()
{
    // Only apply if there are actual padding bits.
    // if(( size_ % IntSize ) == 0 ) {
    //     return;
    // }
    // --> Nope, we have changed the mask to be all-one for its first entry, so that we
    // can avoid the branching here!

    assert( size_ % IntSize < ones_mask_.size() );
    data_.back() &= ones_mask_[ size_ % IntSize ];

    // other versions that might be helpful if i messed up with this little/big endian stuff...
    // first one is slow but definitely works, second one is fast, but might have the same
    // issue as the used version above (which currently works perfectly).
    //~ for (size_t i = size_ % IntSize; i < IntSize; ++i) {
        //~ data_.back() &= ~bit_mask_[i];
    //~ }
    //~ data_.back() &= bit_mask_[size_ % IntSize] - 1;
}

// =============================================================================
//     Dump and Debug
// =============================================================================

std::string Bitvector::dump() const
{
    std::string res = "[" + std::to_string(size_) + "]\n";
    for (size_t i = 0; i < size_; ++i) {
        res += (*this)[i] ? "1" : "0";
        if ((i+1) % 64 == 0) {
            res += "\n";
        } else if ((i+1) % 8 == 0) {
            res += " ";
        }
    }
    return res;
}

std::string Bitvector::dump_int(IntType x) const
{
    std::string res = "";
    for (size_t i = 0; i < IntSize; ++i) {
        res += (x & bit_mask_[i] ? "1" : "0");
        if ((i+1) % 8 == 0) {
            res += " ";
        }
    }
    return res;
}

} // namespace utils
} // namespace genesis

