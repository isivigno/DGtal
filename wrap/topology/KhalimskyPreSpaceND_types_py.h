/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
#ifndef DGTAL_KHALIMSKYPRESPACEND_TYPES_PY_H
#define DGTAL_KHALIMSKYPRESPACEND_TYPES_PY_H

#if defined (_MSC_VER) and !defined(ssize_t)
    // ssize_t is not standard, only posix which is not supported by MSVC
    #define ssize_t ptrdiff_t
#endif

#include "base/Common_types_py.h" // For DGtal::Python::Integer
#include "DGtal/topology/KhalimskyPreSpaceND.h"

namespace DGtal {
    namespace Python {
        using KPreSpace2D = DGtal::KhalimskyPreSpaceND<2, Python::Integer>;
        using PreCell2D = KPreSpace2D::Cell;
        using SPreCell2D = KPreSpace2D::SCell;

        using KPreSpace3D = DGtal::KhalimskyPreSpaceND<3, Python::Integer>;
        using PreCell3D = KPreSpace3D::Cell;
        using SPreCell3D = KPreSpace3D::SCell;
    } // namespace Python
} // namespace DGtal
#endif
