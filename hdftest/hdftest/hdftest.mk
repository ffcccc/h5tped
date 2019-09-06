##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=hdftest
ConfigurationName      :=Debug
WorkspaceConfiguration := $(ConfigurationName)
WorkspacePath          :=/home/fabio/workspace/HUGEF/h5cpp/h5cpp
ProjectPath            :=/home/fabio/workspace/HUGEF/h5cpp/hdftest/hdftest
IntermediateDirectory  :=../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest
OutDir                 :=../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=fabio
Date                   :=04/09/19
CodeLitePath           :=/home/fabio/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=../../h5cpp/build-$(ConfigurationName)/bin/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :=$(IntermediateDirectory)/ObjectsList.txt
PCHCompileFlags        :=
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)../hdfvar $(IncludeSwitch)/usr/include/hdf5/serial 
IncludePCH             := 
RcIncludePath          := 
Libs                   := $(LibrarySwitch)hdfvar $(LibrarySwitch)hdf5_cpp $(LibrarySwitch)hdf5_serial $(LibrarySwitch)hdf5_hl_cpp $(LibrarySwitch)hdf5_serial_hl $(LibrarySwitch)pthread $(LibrarySwitch)z $(LibrarySwitch)sz $(LibrarySwitch)dl $(LibrarySwitch)m $(LibrarySwitch)rt 
ArLibs                 :=  "hdfvar" "hdf5_cpp" "hdf5_serial" "hdf5_hl_cpp" "hdf5_serial_hl" "pthread" "z" "sz" "dl" "m" "rt" 
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)$(WorkspacePath)/../lib $(LibraryPathSwitch)$(WorkspacePath)/build-$(WorkspaceConfiguration)/lib 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(ObjectSuffix) ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(ObjectSuffix) ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: MakeIntermediateDirs $(OutputFile)

$(OutputFile): ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/.d $(Objects) 
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest"
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest"
	@mkdir -p ""../../h5cpp/build-$(ConfigurationName)/bin""

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/.d:
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest"

PreBuild:


##
## Objects
##
../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(ObjectSuffix): ../DataFile.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdftest/DataFile.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_DataFile.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(DependSuffix): ../DataFile.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(DependSuffix) -MM ../DataFile.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(PreprocessSuffix): ../DataFile.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_DataFile.cpp$(PreprocessSuffix) ../DataFile.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(ObjectSuffix): ../main.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdftest/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_main.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(DependSuffix): ../main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(DependSuffix) -MM ../main.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(PreprocessSuffix): ../main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_main.cpp$(PreprocessSuffix) ../main.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(ObjectSuffix): ../hugeDataFile.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdftest/hugeDataFile.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_hugeDataFile.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(DependSuffix): ../hugeDataFile.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(DependSuffix) -MM ../hugeDataFile.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(PreprocessSuffix): ../hugeDataFile.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest/up_hugeDataFile.cpp$(PreprocessSuffix) ../hugeDataFile.cpp


-include ../../h5cpp/build-$(ConfigurationName)/__/hdftest/hdftest//*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r $(IntermediateDirectory)


