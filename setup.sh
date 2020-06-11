#Make sure you have git in your system
git init
git clone https://github.com/mzwang2012/SCauchy
git clone https://github.com/idies/turblib
cd turblib
cp TurbulenceService.h TurbulenceServiceSoap.nsmap soapC.c soapClient.c soapH.h soapStub.h stdsoap2.c stdsoap2.h turblib.c turblib.h ../SCauchy/

:'
TO RUN THE CODE, CHECK README IN SCAUCHY OR VISIT https://github.com/mzwang2012/SCauchy

#-----SOME TIPS------#
TO RUN THE CODE, GO TO SCauchy and run the following commands

-----------------------------------
make all (ignore the warnings)
mkdir checkpoint
./scauchy
----------------------------------

In case you see some errors that you did not see before or do not understand, here is a quick fix tip. Save the checkpoint folders and do the following

rm -r checkpoints
make clean

Then use above mentioned steps and build and run again. Also do this in case you want to run the code again but do not want to start from a checkpoint.

 

'

