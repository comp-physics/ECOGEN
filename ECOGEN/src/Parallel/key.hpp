//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#ifndef INCLUDED_KEY_HPP
#define INCLUDED_KEY_HPP

//! \file      key.hpp
//! \author    B. Dorschner
//! \version   1.0
//! \date      February 15 2019

#include <array>
#include <iomanip>
#include <cmath>
#include "vector.hpp"
#include "bitmasks.hpp"

namespace decomposition
{

template<int Dim=3>
struct Key
{
public: // member types

    template<typename U>
    using vector_type = math::vector<U,Dim>;
    using float_type = double;

    using value_type = bitmask_t::index_t;
    using difference_type = bitmask_t::difference_type;
    using scalar_coordinate_type = bitmask_t::scalar_coordinate_type;
    using coordinate_type = vector_type<scalar_coordinate_type>;;
    using level_type = bitmask_t::level_type;
    using real_scalar_coordinate_type = float_type;
    using real_coordinate_type = vector_type<float_type>;

public: // static

    static constexpr level_type max_level() noexcept
    {
        return bitmask_t::max_level;
    }

    static Key begin(level_type _level) noexcept
    {
        return {bitmask_t::min_arr[_level]};
    }

    static Key end(level_type _level) noexcept
    {
        return {(bitmask_t::min_arr[_level] | bitmask_t::coord_mask)+1u};
    }

    static Key rbegin(level_type _level) noexcept
    {
        return {bitmask_t::max_arr[_level]};
    }

    static Key rend(level_type _level) noexcept
    {
        return {bitmask_t::min_arr[_level]-1u};
    }

    static Key min(level_type _level) noexcept
    {
        return {bitmask_t::min_arr[_level]};
    }

    static Key max(level_type _level) noexcept
    {
        return {bitmask_t::max_arr[_level]};
    }

    static value_type max_index(level_type _level) noexcept
    {
        return (bitmask_t::coord_mask_arr[_level] >> (7+(bitmask_t::max_level-_level)*3));
    }

    static Key at(value_type index, level_type _level) noexcept
    {
        if (index > max_index(_level)) return end(_level);
        index = index << (7+(bitmask_t::max_level-_level)*3);
        return Key{ index | bitmask_t::min_arr[_level] };
    }


    static scalar_coordinate_type max_coordinate(level_type _level) noexcept
    {
        return bitmask_t::max_coord_arr[_level];
    }

   static coordinate_type clamp(const coordinate_type& x, int level) noexcept
   {
    return {{ std::min(std::max(x.x(),0),bitmask_t::max_coord_arr[level]-1),
              std::min(std::max(x.y(),0),bitmask_t::max_coord_arr[level]-1),
              std::min(std::max(x.z(),0),bitmask_t::max_coord_arr[level]-1)}};
   }

   static level_type minimum_level(coordinate_type  x) noexcept
   {
        scalar_coordinate_type m=std::max(std::max(x.x(),x.y()),x.z());
        level_type _level=0;
        while(m>=bitmask_t::max_coord_arr[_level])
        {
            ++_level;
        }
        return _level;
    }

    static bool representable(const coordinate_type& _x, level_type _level)
    {
        scalar_coordinate_type m=std::max(std::max(_x.x(),_x.y()),_x.z());
        scalar_coordinate_type mini=std::min(std::min(_x.x(),_x.y()),_x.z());

        if(m<=bitmask_t::max_coord_arr[_level]-1 && mini>=0)
            return true;

        return false;
    }


public: // Ctors

    Key() noexcept
    : _index(bitmask_t::min_0) {}

    Key(value_type idx) noexcept
    : _index(idx) {}

    Key(coordinate_type x, level_type _level) noexcept
    {
        _index = bitmask_t::lo_mask;
        x = clamp(x,_level);
        value_type _x(x.x()), _y(x.y()), _z(x.z());
        for (auto i=_level; i>0; --i)
        {
            _index |= ((_x & bitmask_t::lo_mask) << ((bitmask_t::max_level-i)*3+7));
            _index |= ((_y & bitmask_t::lo_mask) << ((bitmask_t::max_level-i)*3+8));
            _index |= ((_z & bitmask_t::lo_mask) << ((bitmask_t::max_level-i)*3+9));
            _x >>= 1;
            _y >>= 1;
            _z >>= 1;
        }
        _index |= (static_cast<value_type>(_level) << 2);
    }

    Key(int x, int y, int z, level_type _level) noexcept
    :Key(coordinate_type({x,y,z}), _level) { }

    Key(coordinate_type _x)
    : Key(_x, minimum_level(_x)) { }



    Key(const Key&) = default;
    Key(Key&&) = default;
    Key& operator=(const Key&) & = default;
    Key& operator=(Key&&) & = default;

public: // advance

    Key& operator++() noexcept
    {
        level_type _level = level();
        if (_index >= max(_level)._index) {_index = end(_level)._index; return *this;};
        if (_index <= rend(_level)._index) {_index = min(_level)._index; return *this;};
        _index = (((((bitmask_t::coord_mask_arr[_level] & _index)
                        >> ((bitmask_t::max_level-_level)*3+7)) + 1u)
                 << ((bitmask_t::max_level-_level)*3+7)) | bitmask_t::min_arr[_level]);
        return *this;
    }

    Key operator++(int) noexcept
    {
        Key tmp(*this);
        this->operator++();
        return tmp;
    }

    Key& operator--() noexcept
    {
        level_type _level = level();
        if (_index >= end(_level)._index) {_index = max(_level)._index; return *this;};
        if (_index <= min(_level)._index) {_index = rend(_level)._index; return *this;};
        _index = (((((bitmask_t::coord_mask_arr[_level] & _index)
                        >> ((bitmask_t::max_level-_level)*3+7)) - 1u)
                 << ((bitmask_t::max_level-_level)*3+7)) | bitmask_t::min_arr[_level]);
        return *this;
    }

    Key operator--(int) noexcept
    {
        Key tmp(*this);
        this->operator--();
        return tmp;
    }

    Key& operator +=(int n) noexcept
    {
        level_type _level = level();
        if (n < 0 && _index > (bitmask_t::max_arr[_level]))
        {
            _index = bitmask_t::max_arr[_level];
            ++n;
        }
        else if (_index > (bitmask_t::max_arr[_level]))
            return *this;

        if (n > 0 && _index < (bitmask_t::min_arr[_level]))
        {
            _index = bitmask_t::min_arr[_level];
            --n;
        }
        else if (_index < (bitmask_t::min_arr[_level]))
            return *this;

        if (n==0) return *this;

        difference_type c  = ((bitmask_t::coord_mask_arr[_level] & _index)
                               >> ((bitmask_t::max_level-_level)*3+7));

        if (n < 0)
        {
            if (-n > c) {_index = rend(_level)._index; return *this;};
        }
        else
        {
            const difference_type m  = ((bitmask_t::coord_mask_arr[_level]) >>
                                       ((bitmask_t::max_level-_level)*3+7));
            if (n > m-c) {_index = end(_level)._index; return *this;};
        }

        c+=n;
        _index = ((static_cast<value_type>(c) <<
                 ((bitmask_t::max_level-_level)*3+7)) | bitmask_t::min_arr[_level]);
        return *this;
    }

    Key& operator -=(int n) noexcept
    {
        return this->operator+=(-n);
    }

    Key operator+(int n) const noexcept
    {
        Key tmp(*this);
        return tmp+=n;
    }

    Key operator-(int n) const noexcept
    {
        return this->operator+(-n);
    }

    difference_type operator-(const Key& x) const noexcept
    {
        const auto level_l = level();
        const auto level_r = x.level();
        const auto _level = std::max(level_l,level_r);
        difference_type res(0);
        if (_index >= end(level_l)._index) res += 1;
        else if (_index <= rend(level_l)._index) res -= 1;
        if (x._index >= end(level_r)._index) res -= 1;
        else if (x._index <= rend(level_r)._index) res += 1;
        if (_level == 0) return res;
        const auto l = ((_index) >> (7+3*(bitmask_t::max_level-_level)));
        const auto r = ((x._index) >> (7+3*(bitmask_t::max_level-_level)));
        return res + static_cast<difference_type>(l) - static_cast<difference_type>(r);
    }

public: // queries

    int sib_number(level_type l) const noexcept
    {
        return ((bitmask_t::hi_3_mask >> l*3) & _index) >> (61-l*3);
    }

    level_type level() const noexcept
    {
        return ((bitmask_t::level_mask & _index) >> 2);
    }

    coordinate_type coordinate() const noexcept
    {
        level_type _level = level();
        value_type x(0u), y(0u), z(0u);
        if (_level==0) return {{0,0,0}};
        auto k = (_index >> (7 + 3*(bitmask_t::max_level-_level)));
        for (level_type i=0; i<_level; ++i)
        {
            x |= ((k & bitmask_t::lo_mask) << i);
            k >>= 1;
            y |= ((k & bitmask_t::lo_mask) << i);
            k >>= 1;
            z |= ((k & bitmask_t::lo_mask) << i);
            k >>= 1;
        }
        return {{static_cast<scalar_coordinate_type>(x),
                 static_cast<scalar_coordinate_type>(y),
                 static_cast<scalar_coordinate_type>(z)}};
    }

    real_coordinate_type coordinate_1() const noexcept
    {
        return real_coordinate_type(coordinate())/bitmask_t::max_coord_arr[level()];
    }

    Key neighbor(const coordinate_type& _offset) const noexcept
    {
        const auto c =coordinate();
        const auto l = level();
        const auto cc=_offset+c;
        if(representable(cc,l))
        {
            return Key( cc, level() );
        } else{
            //std::cout<<"non-representable"<<std::endl;
            return end(l);
        }
    }


    Key parent() const noexcept
    {
        const auto _level(level());
        if (_level == 0) return *this;
        return {((bitmask_t::coord_mask_arr[_level-1] & _index) |
                (static_cast<value_type>(_level - 1) << 2))+1u};
    }

    Key level_up_to(level_type Lv ) const noexcept
    {
        if (Lv == 0) return *this;
        return {((bitmask_t::coord_mask_arr[Lv] & _index) |
                (static_cast<value_type>(Lv) << 2))+1u};
    }

    Key child(int i) const noexcept
    {
        const auto _level(level());
        if (_level == bitmask_t::max_level) return *this;
        return {((bitmask_t::coord_mask_arr[_level] & _index) |
                (static_cast<value_type>(i) << ((bitmask_t::max_level-(_level+1))*3+7)) |
                (static_cast<value_type>(_level+1) << 2))+1u};
    }

    int child_number() const noexcept
    {
        const auto _level(level());
        if (_level == 0) return -1;
        return (((bitmask_t::coord_mask_arr[_level] & _index) -
                 (bitmask_t::coord_mask_arr[_level-1] & _index))
                >> (7 + 3*(bitmask_t::max_level-_level)));
    }

    Key equal_coordinate_parent() const noexcept
    {
        Key p(*this);
        if (is_end())
            return p;
        else if (is_rend())
            p = rend(0);
        else
            while (p.child_number()==0)
                p = p.parent();
        return p;
    }

    Key equal_coordinate_child() const noexcept
    {
        Key c(*this);
        if (is_rend())
            return c;
        else if (is_end())
            c = end(bitmask_t::max_level);
        else
            for (auto l=level(); l<bitmask_t::max_level; ++l)
                c = c.child(0);
        return c;
    }

    std::pair<Key,Key> equal_coordinate_range() const noexcept
    {
        return std::make_pair(
            equal_coordinate_parent(),
            (++equal_coordinate_child()).equal_coordinate_parent());
    }

    bool is_end() const noexcept
    {
        return _index >= end(level())._index;
    }


private:

  
    bool is_rend() const noexcept
    {
        return _index <= rend(level())._index;
    }

public: // print

    friend std::ostream& operator<<(std::ostream& os, const Key& t)
    {
        os
        << "{ "
        << "\033[1m";
        for (int i=0; i<bitmask_t::max_level; ++i)
        {
            if (i%2 == 1)
                os << "\033[38;2;159;63;63m";
            else
                os << "\033[38;2;63;127;159m";
            os
            << (((t._index & (bitmask_t::hi_mask >> (i*3)))>0u)?'1':'0')
            << (((t._index & (bitmask_t::hi_mask >> (i*3+1)))>0u)?'1':'0')
            << (((t._index & (bitmask_t::hi_mask >> (i*3+2)))>0u)?'1':'0');
        }
        os << "\033[38;2;127;159;127m";
        for (int i=bitmask_t::max_level*3; i<(bitmask_t::max_level*3+5); ++i)
        {
            os << (((t._index & (bitmask_t::hi_mask >> i))>0u)?'1':'0');
        }
        os
        << "\033[38;2;127;127;127m"
        << (((t._index & (bitmask_t::lo_mask << 1))>0u)?'1':'0')
        << (((t._index & (bitmask_t::lo_mask))>0u)?'1':'0')
        << "\033[0m"
        << " = " << std::setw(22) << std::right << t._index
        << ", level = " << std::setw(2) << std::right << t.level()
        //<< ", coord = " << t.coordinate()
        << ", coord " <<  (t.is_end()?">":(t.is_rend()?"<":"="))
        << " ("  << std::setw(6) << std::right << t.coordinate().x()
        << ", "  << std::setw(6) << std::right << t.coordinate().y()
        << ", "  << std::setw(6) << std::right << t.coordinate().z()
        << ") = " << std::fixed
        << " ("  << std::setw(11) << std::left << std::setprecision(9) << t.coordinate_1().x()
        << ", "  << std::setw(11) << std::left << std::setprecision(9) << t.coordinate_1().y()
        << ", "  << std::setw(11) << std::left << std::setprecision(9) << t.coordinate_1().z()
        << ") }" << std::left; // << std::defaultfloat;
        return os;
    }

    void print_binary(value_type _idx)
    {

        std::cout<< "\033[1m";
        for (int i=0; i<bitmask_t::max_level; ++i)
        {
            if (i%2 == 1)
                std::cout << "\033[38;2;159;63;63m";
            else
                std::cout << "\033[38;2;63;127;159m";
            std::cout
                << (((_idx & (bitmask_t::hi_mask >> (i*3)))>0u)?'1':'0')
                << (((_idx & (bitmask_t::hi_mask >> (i*3+1)))>0u)?'1':'0')
                << (((_idx & (bitmask_t::hi_mask >> (i*3+2)))>0u)?'1':'0');
        }
        std::cout << "\033[38;2;127;159;127m";
        for (int i=bitmask_t::max_level*3; i<(bitmask_t::max_level*3+5); ++i)
        {
            std::cout << (((_idx & (bitmask_t::hi_mask >> i))>0u)?'1':'0');
        }
        std::cout
        << "\033[38;2;127;127;127m"
        << (((_idx & (bitmask_t::lo_mask << 1))>0u)?'1':'0')
        << (((_idx & (bitmask_t::lo_mask))>0u)?'1':'0')
        << "\033[0m";

        std::cout<<std::endl;
    }



public: // members

    value_type _index;
};

// binary operators
template<int Dim>
Key<Dim> operator+(int n, Key<Dim> k) noexcept
{ return k+=n; }

template<int Dim>
Key<Dim> operator-(int n, Key<Dim> k) noexcept
{ return k-=n; }

// relational operators

template<int Dim>
constexpr bool operator==(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index == r._index; }

template<int Dim>
constexpr bool operator!=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index != r._index; }

template<int Dim>
constexpr bool operator<(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index < r._index; }

template<int Dim>
constexpr bool operator<=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index <= r._index; }

template<int Dim>
constexpr bool operator>(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index > r._index; }

template<int Dim>
constexpr bool operator>=(const Key<Dim>& l, const Key<Dim>& r) noexcept
{ return l._index >= r._index; }

// range
template<int Dim>
std::pair<Key<Dim>,Key<Dim>> normalized_range(const Key<Dim>& l, const Key<Dim>& r) noexcept
{
    return std::make_pair(l.equal_coordinate_range().first, r.equal_coordinate_range().second);
}

}

#endif
