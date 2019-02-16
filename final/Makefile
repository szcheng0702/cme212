
# Setup compiler options.
CXX= g++
CXXFLAGS= -std=c++11 -Wall -Wconversion -Wextra
SFMLPATH = /usr/include/SFML
LIBPATH = -L$(SFMLPATH)/lib
INCLUDEPATH = -I$(SFMLPATH)/include
LIBS = -lsfml-graphics -lsfml-window -lsfml-system -lX11 -lGL -lpthread
OBJS = main.o

# Specify the target file, i.e. the executable name.
TARGET= final

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) -o $@ $^ $(LIBPATH) $(LIBS)

# Make pattern rule, to generate an object file from a source code file.
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) -c $<

main.o : main.cpp Viewer.hpp 

.PHONY: clean
clean:
	-rm $(TARGET)
	-rm *.o

#No need to add instructions for building objects here. The SFML library is in 
#the default directory, i.e. /usr/include/, in which g++ compiler will look for
#the needed header files directory. The program already runs without warning 
#under the flags -std=c++11 =Wall =Wconversion =Wextra

