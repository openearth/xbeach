FROM library/centos:7 as builder
RUN yum update -y && \
	yum install epel-release -y 
	
RUN	yum install -y gcc \
	hdf5-devel \
	libtool \
	make \
	netcdf-devel \
	netcdf-fortran \
	netcdf-fortran-devel \
	openmpi-devel \
	python-mako \
	subversion

ENV PATH=/usr/lib64/openmpi/bin:$PATH \
	LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH	

COPY . /opt/xbeach

WORKDIR /opt/xbeach

RUN ./autogen.sh && \
	FCFLAGS=-I/usr/lib64/gfortran/modules ./configure --with-netcdf --without-mpi && \
	make clean && \
	make && \
	make install


FROM library/centos:7

RUN yum update -y && \
	yum install epel-release -y && \
	yum install -y \
	geos \
	git \
	netcdf-fortran-openmpi \
	netcdf-fortran \
	openmpi \
	python34 \
	python34-pip \
	&&  yum clean all

RUN pip3 install mmi

RUN git clone -b mmi https://github.com/openearth/xbeach-mi.git /opt/xbeach-mi
RUN cd /opt/xbeach-mi && python3 setup.py install

ENV PATH=/usr/lib64/openmpi/bin:$PATH \
	LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH \
	PORT=53606

EXPOSE 53606-53620

COPY --from=builder /usr/local/ /usr/local/

WORKDIR /data
CMD mmi-runner --port $PORT --pause xbeach params.txt
