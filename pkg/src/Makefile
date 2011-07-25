
all: rbamtools.so

rbamtools.so: ./rbamtools.o  ./bamtools/BGZF.o ./bamtools/BamReader.o ./bamtools/BamWriter.o
	@echo 'creating rbamtools.so:'
	g++ -fPIC -shared -o"rbamtools.so"  ./rbamtools.o  ./bamtools/BGZF.o ./bamtools/BamReader.o ./bamtools/BamWriter.o -lz

rbamtools.o: rbamtools.cpp
	@echo 'creating rbamtools.o:'
	g++ -I/usr/share/R/include -O3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"rbamtools.d" -MT"rbamtools.d" -o"rbamtools.o" "rbamtools.cpp"

./bamtools/BGZF.o: ./bamtools/BGZF.cpp
	@echo 'Creating BGZF.o:'
	g++ -I/usr/share/R/include -O3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"./bamtools/BGZF.d" -MT"./bamtools/BGZF.d" -o"./bamtools/BGZF.o" "./bamtools/BGZF.cpp"

./bamtools/BamReader.o: ./bamtools/BamReader.cpp
	@echo 'creating BamReader.o:'
	g++ -I/usr/share/R/include -O3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"./bamtools/BamReader.d" -MT"./bamtools/BamReader.d" -o"./bamtools/BamReader.o" "./bamtools/BamReader.cpp"

./bamtools/BamWriter.o: ./bamtools/BamWriter.cpp
	@echo 'creating BamWriter.o:'
	g++ -I/usr/share/R/include -O3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"./bamtools/BamWriter.d" -MT"./bamtools/BamWriter.d" -o"./bamtools/BamWriter.o" "./bamtools/BamWriter.cpp"
clean:
	rm -rf ./rbamtools.o ./rbamtools.d ./bamtools/BamWriter.o ./bamtools/BamWriter.d ./bamtools/BamReader.d ./bamtools/BamReader.o ./bamtools/BGZF.d ./bamtools/BGZF.o
