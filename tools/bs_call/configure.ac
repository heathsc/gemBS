AC_INIT([GEMTools],[1.0])

AC_PROG_CC_STDC
AC_CHECK_LIB([z], [gzopen], [HAVE_ZLIB=1])
AC_SUBST([HAVE_ZLIB])
AC_CHECK_LIB([bz2], [BZ2_bzReadOpen], [HAVE_BZLIB=1])
AC_SUBST([HAVE_BZLIB])

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

AC_CONFIG_FILES([gt/Makefile.mk])
AC_CONFIG_FILES([src/Makefile.mk])
AC_OUTPUT
