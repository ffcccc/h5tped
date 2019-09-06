##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=hdfvar
ConfigurationName      :=Debug
WorkspaceConfiguration := $(ConfigurationName)
WorkspacePath          :=/home/fabio/workspace/HUGEF/h5cpp/h5cpp
ProjectPath            :=/home/fabio/workspace/HUGEF/h5cpp/hdfvar/hdfvar
IntermediateDirectory  :=../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar
OutDir                 :=../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar
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
OutputFile             :=../../h5cpp/build-$(ConfigurationName)/lib/lib$(ProjectName).a
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :=$(IntermediateDirectory)/ObjectsList.txt
PCHCompileFlags        :=
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). $(IncludeSwitch)/usr/include/hdf5/serial 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g $(Preprocessors)
CFLAGS   :=  -g $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(ObjectSuffix) ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(ObjectSuffix) ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: MakeIntermediateDirs ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/$(OutputFile)

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/$(OutputFile): $(Objects)
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar"
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(AR) $(ArchiveOutputSwitch)$(OutputFile) @$(ObjectsFileList) $(ArLibs)
	@echo rebuilt > $(IntermediateDirectory)/hdfvar.relink

MakeIntermediateDirs:
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar"
	@mkdir -p ""../../h5cpp/build-$(ConfigurationName)/lib""

./$(WorkspaceConfiguration):
	@mkdir -p "../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar"

PreBuild:


##
## Objects
##
../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(ObjectSuffix): ../hdfTPed.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdfvar/hdfTPed.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_hdfTPed.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(DependSuffix): ../hdfTPed.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(DependSuffix) -MM ../hdfTPed.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(PreprocessSuffix): ../hdfTPed.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfTPed.cpp$(PreprocessSuffix) ../hdfTPed.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(ObjectSuffix): ../hdfvar_I.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdfvar/hdfvar_I.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_hdfvar_I.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(DependSuffix): ../hdfvar_I.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(DependSuffix) -MM ../hdfvar_I.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(PreprocessSuffix): ../hdfvar_I.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_hdfvar_I.cpp$(PreprocessSuffix) ../hdfvar_I.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(ObjectSuffix): ../var_utils.cpp ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/fabio/workspace/HUGEF/h5cpp/hdfvar/var_utils.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/up_var_utils.cpp$(ObjectSuffix) $(IncludePath)
../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(DependSuffix): ../var_utils.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(ObjectSuffix) -MF../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(DependSuffix) -MM ../var_utils.cpp

../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(PreprocessSuffix): ../var_utils.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar/up_var_utils.cpp$(PreprocessSuffix) ../var_utils.cpp


-include ../../h5cpp/build-$(ConfigurationName)/__/hdfvar/hdfvar//*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r $(IntermediateDirectory)


