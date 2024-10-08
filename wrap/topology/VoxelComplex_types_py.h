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

#ifndef DGTAL_VOXELCOMPLEX_TYPES_PY_H
#define DGTAL_VOXELCOMPLEX_TYPES_PY_H

#if defined (_MSC_VER) and !defined(ssize_t)
    // ssize_t is not standard, only posix which is not supported by MSVC
    #define ssize_t ptrdiff_t
#endif

#include "topology/KhalimskySpaceND_types_py.h" // For KSpace3D
#include "topology/CubicalComplex_types_py.h" // For CellMap3D
#include "DGtal/topology/VoxelComplex.h"

namespace DGtal {
    namespace Python {
        using VoxelComplex = DGtal::VoxelComplex<Python::KSpace3D, Python::CellMap3D>;
    } // namespace Python
} // namespace DGtal
#endif
