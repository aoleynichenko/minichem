#input file for autoconf


#0.8 is the package version, it should be changed to correct value
AC_INIT([minichem],[0.8],[alexvoleynichenko@gmail.com])


#is it portable??...
AC_CONFIG_AUX_DIR([/usr/share/libtool/build-aux])

AC_CANONICAL_HOST


#support only x86 and x86-64 - compatible CPUs (Intel and AMD)
case "${build}" in
	i[[3456]]86*)
		CXXFLAGS+="-march=i386 " #is it good? march may be i[456]86
		;;
	x86_64*)
		CXXFLAGS+="-march=x86-64 "
		;;
	*)
		AC_MSG_ERROR([Unsupported arch $build])
esac


case "${host_os}" in 
	linux* )
		CXXFLAGS+="-D__linux__ -D__UNIX_LIKE_OS "
		;;
	*BSD*)
		CXXFLAGS+="-D__linux__ -D__UNIX_LIKE_OS "
		;;
	darwin* )
		CXXFLAGS+="-D__APPLE__ -D__UNIX_LIKE_OS "
		;;
	cygwin | windows* | mingw32*) 
		;;
	* )
		AC_MSG_ERROR([Unsupported host OS $host_os])
		;;
esac

AC_OUTPUT
[
sed "s/.*export CXXFLAGS.*/CXXFLAGS+=$CXXFLAGS\n&/" -i Makefile
]
