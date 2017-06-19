dnl
dnl				acx_urbi_resolve_dir.m4
dnl
dnl Enable a reasonable set of optimization flags for the C compiler. 
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyrights:		(C) 2007-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
dnl				(C) 2007 Akim Demaille (original version)
dnl
dnl	License:		GPL
dnl
dnl	AstrOmatic software is free software: you can redistribute it and/or
dnl	modify it under the terms of the GNU General Public License as
dnl	published by the Free Software Foundation, either version 3 of the
dnl	License, or (at your option) any later version.
dnl	AstrOmatic software is distributed in the hope that it will be useful,
dnl	but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl	GNU General Public License for more details.
dnl	You should have received a copy of the GNU General Public License
dnl	along with AstrOmatic software.
dnl	If not, see <http://www.gnu.org/licenses/>.
dnl
dnl	Last modified:		27/12/2011
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_URBI_RESOLVE_DIR
dnl
dnl Return a directory with all inner variables expanded.
dnl Based on a macro kindly provided by Akim Demaille  <akim@lrde.epita.fr>.

# URBI_RESOLVE_DIR_PREPARE
# ------------------------
# Define urbi_resolve_dir.
m4_defun([URBI_RESOLVE_DIR_PREPARE],
[# PATH urbi_resolve_dir(DIR)
# --------------------------
# 
urbi_resolve_dir ()
{
  ac_$0_dir=$[]1
  ac_$0_res=
  ac_$0_prefix_NONE=
  ac_$0_exec_prefix_NONE=
  test "x$prefix" = xNONE &&
     ac_$0_exec_prefix_NONE=yes &&
     prefix=$ac_default_prefix
  test "x$exec_prefix" = xNONE &&
     ac_$0_exec_prefix_NONE=yes &&
     exec_prefix=$prefix
  while true
  do
    eval ac_$0_res="$ac_$0_dir"
    if test x"$ac_$0_dir" = x"$ac_$0_res"; then
      break
    fi
    ac_$0_dir=$ac_$0_res
  done
  test "$ac_$0_prefix_NONE" && prefix=NONE
  test "$ac_$0_exec_prefix_NONE" && exec_prefix=NONE
  echo "$ac_$0_res"
}
])


# PATH URBI_RESOLVE_DIR(DIR)
# --------------------------
# Return the DIR with all inner variables expanded.
AC_DEFUN([URBI_RESOLVE_DIR],
[AC_REQUIRE([URBI_RESOLVE_DIR_PREPARE])dnl
urbi_resolve_dir '$1'[]dnl
])


## Local Variables:
## mode: autoconf
## End:
