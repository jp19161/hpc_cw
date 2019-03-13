##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=HPC_Coursework
ConfigurationName      :=Debug
WorkspacePath          :=/home/jp1916/HPC_Coursework
ProjectPath            :=/home/jp1916/HPC_Coursework/HPC_Coursework
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Jinay Patel
Date                   :=13/03/19
CodeLitePath           :=/home/jp1916/.codelite
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
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="HPC_Coursework.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). $(LibraryPathSwitch)/home/jp1916/Downloads/blas-3.8.0.tgz 

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
Objects0=$(IntermediateDirectory)/matlabcode.cpp$(ObjectSuffix) $(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) $(IntermediateDirectory)/Model.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/matlabcode.cpp$(ObjectSuffix): matlabcode.cpp $(IntermediateDirectory)/matlabcode.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jp1916/HPC_Coursework/HPC_Coursework/matlabcode.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/matlabcode.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/matlabcode.cpp$(DependSuffix): matlabcode.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/matlabcode.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/matlabcode.cpp$(DependSuffix) -MM matlabcode.cpp

$(IntermediateDirectory)/matlabcode.cpp$(PreprocessSuffix): matlabcode.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/matlabcode.cpp$(PreprocessSuffix) matlabcode.cpp

$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix): Burgers.cpp $(IntermediateDirectory)/Burgers.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jp1916/HPC_Coursework/HPC_Coursework/Burgers.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Burgers.cpp$(DependSuffix): Burgers.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Burgers.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Burgers.cpp$(DependSuffix) -MM Burgers.cpp

$(IntermediateDirectory)/Burgers.cpp$(PreprocessSuffix): Burgers.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Burgers.cpp$(PreprocessSuffix) Burgers.cpp

$(IntermediateDirectory)/Model.cpp$(ObjectSuffix): Model.cpp $(IntermediateDirectory)/Model.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/jp1916/HPC_Coursework/HPC_Coursework/Model.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Model.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Model.cpp$(DependSuffix): Model.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/Model.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/Model.cpp$(DependSuffix) -MM Model.cpp

$(IntermediateDirectory)/Model.cpp$(PreprocessSuffix): Model.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Model.cpp$(PreprocessSuffix) Model.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


