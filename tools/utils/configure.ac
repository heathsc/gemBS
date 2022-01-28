AC_INIT([gemBS_utils],[1.0])

AC_PROG_CC_STDC
AC_ARG_WITH(htslib,
[AS_HELP_STRING([--with-htslib=PATH],
		[specify prefix directory for installed htslib library.])])
if test "x$with_htslib" != "x"; then
      HTSLIBS="$with_htslib/libhts.a -lcurl -lcrypto -llzma"
		HTSINC="-I$with_htslib"
else
      HTSLIBS="-lhts -lcurl -lcrypto -llzma"
		HTSINC=""
fi

AC_SUBST([HTSLIBS])
AC_SUBST([HTSINC])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
