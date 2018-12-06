/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * for internal use only; used to switch between the different matrix libraries
 */

#ifdef STATS_WITH_MATRIX_LIB
namespace mat_ops {
    
    #include "mat_ops/get_mem_ptr.hpp"
    #include "mat_ops/n_cols.hpp"
    #include "mat_ops/n_elem.hpp"
    #include "mat_ops/n_rows.hpp"

    #include "mat_ops/accu.hpp"
    #include "mat_ops/chol.hpp"
    #include "mat_ops/cumsum.hpp"
    #include "mat_ops/det.hpp"
    #include "mat_ops/eye.hpp"
    #include "mat_ops/fill.hpp"
    #include "mat_ops/get_row.hpp"
    #include "mat_ops/inv.hpp"
    #include "mat_ops/log_det.hpp"
    #include "mat_ops/mean.hpp"
    #include "mat_ops/repmat.hpp"
    #include "mat_ops/solve.hpp"
    #include "mat_ops/spacing.hpp"
    #include "mat_ops/trace.hpp"
    #include "mat_ops/trans.hpp"
    #include "mat_ops/var.hpp"
    #include "mat_ops/zeros.hpp"

}
#endif
