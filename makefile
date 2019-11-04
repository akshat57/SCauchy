CC     = gcc -g
FC     = gfortran -g
RM     = rm -f
CFLAGS = -Wall -g -I..
CGLIBS = mtprng.o stdtypes.o
LDLIBS = 

OBJ =	soapC.o \
	soapClient.o \
	stdsoap2.o \
        turblib.o

all: strack

strack : $(OBJ) strack.o  $(CGLIBS)
	 $(FC) -o $@ $(OBJ) strack.o $(CGLIBS) $(LDLIBS) 


strack.o : strack.f90 mtprng.o stdtypes.o 
	$(FC) -c strack.f90

mtprng.o : mtprng.f90 stdtypes.o
	$(FC) -c mtprng.f90

stdtypes.o : stdtypes.f90
	$(FC) -c stdtypes.f90

stdsoap2.o: stdsoap2.c
	$(CC) $(CFLAGS) -c $<

static_lib: $(OBJ)
	ar rcs libJHTDB.a $(OBJ)

install: static_lib
	$(MKDIR) $(JHTDB_PREFIX)/include
	$(MKDIR) $(JHTDB_PREFIX)/lib
	$(CP) *.h $(JHTDB_PREFIX)/include/
	$(CP) libJHTDB.a $(JHTDB_PREFIX)/lib/

# Regenerate the gSOAP interfaces if required
TurbulenceService.h : wsdl

# Update the WSDL and gSOAP interfaces
wsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

testwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://test.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhdtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhdtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

devwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://dev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

mhddevwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://mhddev.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

prodtestwsdl:
	wsdl2h -o TurbulenceService.h -n turb -c "http://prodtest.turbulence.pha.jhu.edu/service/turbulence.asmx?WSDL" -s
	soapcpp2 -CLcx -2 -I.:$(SOAP_INCLUDE_DIR) TurbulenceService.h

clean:
	$(RM) *.o *.exe turbf turbc mhdc mhdf channelc channelf mixingc mixingf compiler_flags

spotless: clean
	$(RM) soapClient.c TurbulenceServiceSoap.nsmap soapH.h TurbulenceServiceSoap12.nsmap soapStub.h soapC.c TurbulenceService.h

.SUFFIXES: .o .c .x

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: force
compiler_flags: force
	echo '$(CFLAGS)' | cmp -s - $@ || echo '$(CFLAGS)' > $@

$(OBJ): compiler_flags

