# These xxx markers are to be replaced by git_build_rpm
Name:           xxx
Version:        xxx

Release:        1%{?dist}
Summary:        Tiny linear algebra library for small matrices

License:        LGPL
URL:            https://github.jpl.nasa.gov/maritime-robotics/libminimath/
Source0:        https://github.jpl.nasa.gov/maritime-robotics/libminimath/archive/%{version}.tar.gz#/%{name}-%{version}.tar.gz

BuildRequires:  mrbuild >= 0.40

%description
This exists because for tiny matrices the overheads of BLAS and LAPACK are
significant

%package devel
Summary: Tiny linear algebra library for small matrices
%description devel
This exists because for tiny matrices the overheads of BLAS and LAPACK are
significant

%prep
%setup -q

%build
make %{?_smp_mflags}

%install
rm -rf $RPM_BUILD_ROOT
%make_install

%post -p /sbin/ldconfig
%postun -p /sbin/ldconfig

%files devel
%doc README
%{_includedir}/*
