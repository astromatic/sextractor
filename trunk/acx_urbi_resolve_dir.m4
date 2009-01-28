

# URBI_RESOLVE_DIR_PREPARE
# ------------------------
# Define urbi_resolve_dir.
m4_defun([URBI_RESOLVE_DIR_PREPARE],
[# PATH urbi_resolve_dir(DIR)
# --------------------------
# Return the DIR with all inner variables expanded.
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
    if test x"$ac_$0_dir" == x"$ac_$0_res"; then
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
