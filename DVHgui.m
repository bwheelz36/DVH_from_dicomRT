function varargout = DVHgui(varargin)
%written by Brendan Whelan 2011.
%%
%This code draws heavily on code written by will ansbacher in 2005 for the 
%epidose program. It also utilises a function written by Karl Otto.

%GUI to create DVHs from input dicom-rt files. Input files must be a dose
%file and a structure file. 

%FUTURE WORK:
% Incorporate the ability to select multiple dose files and loop over them
% Incorporate the option to compute real volumes for the DVHs
% add status bars
%%
% Last Modified by GUIDE v2.5 29-Oct-2013 12:43:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DVHgui_OpeningFcn, ...
                   'gui_OutputFcn',  @DVHgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DVHgui is made visible.
function DVHgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DVHgui (see VARARGIN)
handles.output = hObject;

%set initial states
savePath=pwd;
set(handles.edit4,'String',savePath);
handles.DVHtype='cum';
% %set up the 'cum' and 'dif' button selection function
set(handles.uipanel2,'SelectionChangeFcn',@DVHtype_buttongroup_SelectionChangeFcn);
% % Update handles structure
 guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.

function varargout = DVHgui_OutputFcn(hObject, eventdata, handles) 
%%
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function listbox1_Callback(hObject, eventdata, handles,contours)
%%
index_selected = get(hObject,'Value');
handles.index_selected=index_selected;
guidata(hObject, handles);


function listbox1_CreateFcn(hObject, eventdata, handles)
%%
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit1_Callback(hObject, eventdata, handles)
%%
disp('note that default save path can be changed at line 59 &/or 149')
input = get(hObject,'String');

if isempty(input)
   savePath=pwd;
   set(handles.edit1,'String',savePath);
end
guidata(hObject, handles);


function edit1_CreateFcn(hObject, eventdata, handles)
%%

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton1_Callback(hObject, eventdata, handles)

if isfield(handles,'dosefiles_selected')
    numDoseFiles=numel(handles.dosefiles_selected);
    DoseFiles=handles.doseFiles(handles.dosefiles_selected);
    dicompath=handles.dicompath;
else
    
    DoseFiles=handles.doseFiles;
    numDoseFiles=numel(handles.doseFiles);
    dicompath=handles.dicompath;
end

for jj=1:numDoseFiles
    
        Edosefile=DoseFiles{jj};
  
% get the dose data (w. ansbacher)
DoseHeader = dicominfo(fullfile(dicompath,Edosefile));
DoseID=DoseHeader.PatientID;
Dose    = dicomread(DoseHeader);
Dose   = double(squeeze(Dose))*(DoseHeader.DoseGridScaling);
if isempty(Dose)
    error('This dicom dose file doesnt seem to have any dose in it. I have no idea why why this happens but sometimes it does....');
end
%%
%Build a coordinate system (w. ansbacher)
resXY  = DoseHeader.PixelSpacing;          % x and y resolution in mm
XYZmin = DoseHeader.ImagePositionPatient;  % x, y, z start of image
X = XYZmin(1) + resXY(1)*[ 0 : double(DoseHeader.Width  - 1) ];  % vector in horiz dir
Y = XYZmin(2) + resXY(2)*[ 0 : double(DoseHeader.Height - 1) ];  % vector in vert dir
Z = XYZmin(3) + DoseHeader.GridFrameOffsetVector;
Coords{1}=X;    Coords{2}=Y;    Coords{3}=Z;
 handles.Coords=Coords;
 handles.Dose=Dose;
 handles.DoseID=DoseID;
  
%update handles structure
guidata(hObject, handles);
%% 
%Interpolate the contours and the dose onto the same coords. This was
%written by Will Ansbacher
Dose=handles.Dose;
STRUCT=handles.STRUCT;
cont_file=handles.cont_file;
if ~isequal(handles.structureID,handles.DoseID)
    warning('patient ID from dose file and contour file does not match')
   
end
Coords=handles.Coords;

[Dose Z]=InterpolateDoseContours(STRUCT,Dose,cont_file,Coords);

handles.Coords{3}=Z;
handles.Dose=Dose;


%%
for i=1:numel(handles.index_selected)   

    savePath=strcat(get(handles.edit4,'String'),'\',DoseFiles{jj});
   
if ~isdir(savePath)
    mkdir(savePath)
end
contour=handles.STRUCT{handles.index_selected(i)};
DVHgenerate(contour,handles.Dose,handles.Coords,savePath,handles.DVHtype)
end
fprintf('\n finished \n');
end


function DVHtype_buttongroup_SelectionChangeFcn(hObject, eventdata)
%%
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton3'
        handles.DVHtype='cum';   
    case 'radiobutton4'
        handles.DVHtype='dif';
    otherwise
        %(don't actually think this line is needed but better safe...)
        handles.DVHtype='cum';    
end
%updates the handles structure
guidata(hObject, handles);


function pushbutton2_Callback(hObject, eventdata, handles)
temp='select all dose files';
[Edosefile, dicompath]=uigetfile('V:\Brendan\ED study\DICOM data\DICOM export\IMRT plans\NEW PLANS!\*.dcm',temp,'MultiSelect','on');
if ischar(class(Edosefile)) %to deal with the annoying fact that uigetfile returns characters instead of cells if one item selected
    Edosefile=cellstr(Edosefile);
end
 %Populate ListBox of structures
 
 set(handles.listbox2,'String',Edosefile);
 handles.doseFiles=Edosefile;
 handles.dicompath=dicompath;
 guidata(hObject, handles);
 %later will have to put in a routine for the situation where not all
 %structures are the same... 
 

% %%


function pushbutton3_Callback(hObject, eventdata, handles)
%%
%get structure file
temp='select structure file';
[cont_file, dicompath]=uigetfile('V:\Brendan\ED study\DICOM data\DICOM export\IMRT plans\NEW PLANS!\*.dcm',temp);

%%
%GET CONOTURS FROM RTSS FILE
%written by Karl Otto
fprintf('\n getting the contour data... \n');
structure_info=dicominfo(fullfile(dicompath,cont_file)); 
structureID=structure_info.PatientID;
orientation = 'HFS'; % (not relevant as I search for each contour anyway)
contShift=[0,0,0]'; %this is an offset; not relevant for me because images have exactly the same orientation
[STRUCT,patient_data] =convert_structures_epi(structure_info,orientation,contShift);
%update current dose fiel box
handles.STRUCT=STRUCT;
handles.structureID=structureID;
handles.cont_file=cont_file;
set(handles.edit3,'String',cont_file);
%%
%get a list of contour names to put in the list box

  for j=1:numel(STRUCT)
      contours{j,1}=(STRUCT{1,j}.name);
  end
set(handles.listbox1,'String',contours);
% update hadnles structure
guidata(hObject, handles);

function edit2_Callback(hObject, eventdata, handles)
%%
input = get(hObject,'String');

guidata(hObject, handles);

function edit2_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
%%
input = get(hObject,'String');

guidata(hObject, handles);

function edit3_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function listbox2_Callback(hObject, eventdata, handles)
handles.dosefiles_selected = get(hObject,'Value');
guidata(hObject, handles);


function listbox2_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton5_Callback(hObject, eventdata, handles)
handles.SaveFolder=uigetdir;
set(handles.edit4,'String',handles.SaveFolder);
guidata(hObject, handles);

function edit4_Callback(hObject, eventdata, handles)
handles.SaveFolder = get(hObject,'String');
guidata(hObject, handles);

function edit4_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [Edose Z]=InterpolateDoseContours(STRUCT,Edose,cont_file,Coords)
X=Coords{1};    Y=Coords{2};    Z=Coords{3};
%%%%%%%%%%            Author:  W. Ansbacher           %%%%%%%%%%
%%%%%%%%%%   Copyright BC Cancer Agency Branch 2005   %%%%%%%%%%
 fprintf('\n Interpolating contours...');
  %find the z-plane spacing for the contours in the structures
 x = 1;
 y = 0;
    while y == 0
        try
            ZisInCont= STRUCT{x}.contour{1}(3,1); % first Z in the first structure
            n = 2;
            ZcontIncr= abs(STRUCT{x}.contour{n}(3,1) - STRUCT{x}.contour{1}(3,1));
                while ZcontIncr == 0    %if 1st two contours are on the same plane
                    n = n+1;
                    ZcontIncr= abs(STRUCT{x}.contour{n}(3,1) - STRUCT{x}.contour{1}(3,1));
                end
            y = 1;
        catch
            x = x+1;
        end
        if x > size(STRUCT,2)
            fprintf('There was an error calculating the z slice thickness. The script has aborted./n')
            diary off
            break
        end
    end
     fprintf('\n Patient contour file  : %s',cont_file);
    Znew1= ZisInCont + ZcontIncr* ceil(( Z(1) - ZisInCont)/ZcontIncr);
    Znew2= ZisInCont + ZcontIncr*floor(( Z(end) - ZisInCont)/ZcontIncr);
    %interpolate the dose matrix on to this new base
    %% Comment out old code (2 lines below)
    %tempE = Edose; temp = [Znew1:ZcontIncr:Znew2]'; %temp is a col vect like Z
    %Edose =interp3(X,Y,Z, tempE ,X,Y,temp);
    
    % What is changed: According to Google, the input for the function
    % interp3 needs proper meshgrid() or ndgrid() data, instead of simple
    % column vectors as used before.
    
    X1 = X; Y1 = Y; % Copy values of original column vectors to X1, Y1. These will be used as parameters for meshgrid() 
    
    tempE = Edose; temp = [Znew1:ZcontIncr:Znew2]'; %temp is a col vect like Z
    if isempty(temp)
       temp = [Znew1:-ZcontIncr:Znew2]'; % dodgy way to account for reverse orientation
    end
        
    
    % Generate meshgrid data => X,Y,Z will each be 3D vectors instead of
    % column vectors
    [X,Y,Z] = meshgrid(X1,Y1,Z);
    [X2,Y2,temp] = meshgrid(X1,Y1,temp);
    
    % Call interp3 with these 3D vectors
    if isequal(X,X2) && isequal(Y,Y2) && isequal(Z,temp)
        % no interpoliation needed, and apparently in newer versions of
        % matlab this causes crahses, go matlab you idiot.
        Edose = tempE;
    else
        Edose =interp3(X,Y,Z, tempE ,X2,Y2,temp);
    end
    
    % reshape temp to a column vector and pass values to Z. The column-vector structure
    % is needed to be consistent with the rest of the code
    

    Z = reshape(temp(1,1,:),1,size(temp(1,1,:),3));
    
    
function [STRUCT,patient_data] =convert_structures_epi(dicom_data,orientation, xyz_shift);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creator and owner: Karl Otto
% Copyright 2004
%modified 20081114 W.Ansbacher for Epidose input needs
% xyz_shift is the 3-ROW  x,y,z offset between verification and patient dicom origins
% NB: untried in the 'HFP' orientation!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dicom_struct=struct2cell(dicom_data.ROIContourSequence);
dicom_struct_names=struct2cell(dicom_data.StructureSetROISequence);
num_structures=size(dicom_struct,1);
patient_data.patient_name=dicom_data.PatientName.FamilyName;
patient_data.patient_number=dicom_data.PatientID;
for i=1:num_structures
    STRUCT{i}.name=dicom_struct_names{i}.ROIName;

    if isfield(dicom_struct{i},'ContourSequence')
        dicom_contours=struct2cell(dicom_struct{i}.ContourSequence);
        num_contours=size(dicom_contours,1);

        %if isequal(orientation,'HFS')  %Head First Supine
            for j=1:num_contours
                % extend xyz_shift(3x1) to the # of columns in the contour
                shiftMatrix = repmat(xyz_shift, 1, dicom_contours{j}.NumberOfContourPoints);
                % reshape contour to this format and add the offset
                STRUCT{i}.contour{j} = reshape(dicom_contours{j}.ContourData, size(shiftMatrix)) ...
                    + shiftMatrix;
                
            end

    end
end

function DVHgenerate(contour,Dose,Coords,savePath,option)
%Brendan Whelan 2011
%function to generate a contour based a Dose grid and contour grid which
%have been interpolated onto the same coordinate system. The code in this
%function has been heavily influenced by code written by w. ansbacher
X=Coords{1};    Y=Coords{2};    Z=Coords{3};

cont_name=contour.name;
contours=contour.contour;
fprintf('\n generating DVH for %s \n',cont_name);
ROI=false(size(Dose)); % this starts off as bunch of zeros the size of dose
%compute voxel volume
%Xdim=abs(X(2)-X(1)),    Ydim=abs(Y(2)-Y(1)),    Zdim=abs(Z(2)-Z(1))
%VoxelVolume=Xdim*Ydim*Zdim 
% we want to fill it in with ones inside the contour
    
    for j=1:length(contours)
        
        Zcoord=contours{1,j}(3,1); %this is where the contour should go on the dose file
        epsilon=1e-3; % this is a tolerance which takes into account any difference between contour location and dose location. It should be made as small as possible.
        dum=(abs(Z-Zcoord)<epsilon);
        use=find(dum);
        
        if isempty(use) 
            %Either contour isn't in dose region, or contours and dose do
            %not line up within tolerance. First check to see if contour is
            %outside dose region.
                ZgridSpacing=Z(2)-Z(1);
                 if ~round(mod(abs(Z(1)-Zcoord),ZgridSpacing))
                     warning('contours for this structure appear to extend beyond dose calc region. Extending dose calc region with the assumption of zero dose')
                     %this isn't the nicest or fastest way to code this but I don't
                     %really want to encourage people to calcualte DVHs for
                     %structures outside dose region so will keep it this
                     %way
                     
%                      Znew=[min(Z(:)): ZgridSpacing: Zcoord]';
%                      DoseExtend=size(Znew);
                     Dose(:,:,size(Z(:)) + 1)=0;
                     use=(numel(Z(:))+1); %it doesn't matter where the contour should actually go as Dose=0 anyway                 
                 
                 else
                    dbstop in DVHgui at 436 % this is intended to put a break point at the next statment, obviously any changes in position will stop it working
                    error('contour could not be located in dose calc region, entering debug mode')
                  
                    
                    %if you've arrived here, it means that the code has not
                    %found the contour location Zcoord anywhere within the
                    %dose calc region Z. I would suggest you go back
                    %through the code and try and work out why this is
                    %happening.
                 end
        end
        %Create ROI
        x =contours{1,j}(1,:);
        y = contours{1,j}(2,:);
        
        %CREATE ROI MATRIX (ones in ROI zeros elsewhere)
        BW = roipoly(X,Y,Dose(:,:,1),x,y); %Dose is a dummy image
       % diamond = [0,1,0;1,1,1;0,1,0]; %dilation matrix
       % BW = imdilate(BW, diamond); %this dilates the dose file inwards
       % one pixel. However, I don't think it is needed.
       
        ROI(:,:,use)=BW;        
    end
            selected=find(ROI);
            maxdose = max(Dose(selected));
            diffHist = histc(Dose(selected), [0:.05:maxdose]);
     %%
    switch option
        case 'dif'
             Hist=diffHist;    
        case 'cum'
            Hist = 100*(1 - cumsum(diffHist)/numel(Dose(selected)));    
                                  
    end
            name=char(cont_name);    filename=fullfile(savePath,name); 
            A(1,:)=0:.05:maxdose;
            A(2,:)=Hist;
            A=A';
            try
            xlswrite(filename,A);
            catch err
                display('likely issue is the naming of the structures...')
                rethrow(err)
                
            end
    

    
    
