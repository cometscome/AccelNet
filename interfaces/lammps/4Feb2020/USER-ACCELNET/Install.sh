#!/bin/sh

# Install/uninstall this package in a LAMMPS 4Feb20 source tree.
# mode = 0/1/2 for uninstall/install/update

mode=$1
LC_ALL=C
export LC_ALL

action () {
  if test "$mode" = 0; then
    rm -f "../$1"
  elif ! cmp -s "$1" "../$1"; then
    cp "$1" ..
    if test "$mode" = 2; then
      echo "  updating src/$1"
    fi
  fi
}

for file in *.cpp *.h; do
  test -f "$file" && action "$file"
done

if test ! -e ../Makefile.package; then
  exit 0
fi

gfortran_library=$(gfortran -print-file-name=libgfortran.a 2>/dev/null)
if test -n "$gfortran_library" && test "$gfortran_library" != libgfortran.a; then
  gfortran_library_dir=$(dirname "$gfortran_library")
else
  gfortran_library_dir=
fi

# Keep repeated package install/update operations idempotent.
if test -n "$gfortran_library_dir"; then
  sed -i.bak -e "s|-L$gfortran_library_dir ||g" ../Makefile.package
fi

if test "$mode" = 1; then
  sed -i.bak -e 's/-laccelnet -lAccelNetDescriptors -lgfortran -lquadmath //' ../Makefile.package
  sed -i.bak -e 's/[^ \t]*accelnet[^ \t]* //' ../Makefile.package
  sed -i.bak -e 's|^PKG_INC =[ \t]*|&-I../../lib/accelnet/include |' ../Makefile.package
  sed -i.bak -e "s|^PKG_PATH =[ \t]*|&-L../../lib/accelnet/lib ${gfortran_library_dir:+-L$gfortran_library_dir }|" ../Makefile.package
  sed -i.bak -e 's|^PKG_LIB =[ \t]*|&-laccelnet -lAccelNetDescriptors -lgfortran -lquadmath |' ../Makefile.package
elif test "$mode" = 0; then
  sed -i.bak -e 's/-laccelnet -lAccelNetDescriptors -lgfortran -lquadmath //' ../Makefile.package
  sed -i.bak -e 's/[^ \t]*accelnet[^ \t]* //' ../Makefile.package
fi

rm -f ../Makefile.package.bak
