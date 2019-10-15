classdef EClab_CV
    % this class define some methoods for using CV 
    % version 1.1 latest update 10/10/19
    properties
        dataTable
        E
        I
        N_cycles
        speed
        ImaxOx
        EmaxOx
        ImaxRed
        EmaxRed
        peaks
    end
    
    methods (Static)
        
        function [ dataWithoutOutliers ] = rmoutliers( data )
            % remove outliers in data
                threshold = 3 * std( data );
                validRange = mean( data ) + [-1 1] * threshold;
                dataWithoutOutliers = data( data >= validRange(1) & data <= validRange(2) );

        end
        
        
        function [varargout]=plot_(a,pairs)
             % enable to plot an arrey of DPVs 
             % if output defined return plot handle
            
            if nargout==1
                eval(['varargout{1}=plot(a.E,a.I' pairs ')'])
            else
                eval(['plot(a.E,a.I' pairs ')'])
            end
            xlabel('E [V]')
            ylabel('I [A]')
        end
        
        
        
        function r = minus_(a,b)
%             % enable to do difrrentiation element wise
            if length(a.E)==length(b.E) && sum(abs(a.E-b.E)<(10^-3))==length(a.E)
                r = a;
                r.I = a.I - b.I;
                r = getExtrimom(r);
            else 
%                 % this part interpulate for not even spaced samples
                if abs(max(a.E)-max(b.E))<10^-3 && abs(min(a.E)-min(b.E))<10^-3
                    temp=evenup([a,b]);   
                    a=temp(1);
                    b=temp(2);
                  
                    r = a;
                    r.I = a.I - b.I;
                    r = getExtrimom(r);
                else
                    error('E range must be the same')
                end
            end
        end
        
       
        
        function a = CV_alloc(varargin)
            % allow CV pre allocation
            warning('off','all')
            dim_str='';
            for i=1:nargin
                if ~(i==nargin)
                    dim_str=[dim_str num2str(varargin{i}) ','];
                else
                    dim_str=[dim_str num2str(varargin{i})];
                end
            end
            eval(['a(' dim_str ') = EClab_CV();']);
            warning('on','all')
        end
        
    end
    methods
        
        function obj = EClab_CV(fileName)
            % initiation function upon calling to the class
            % this part get the properties from the EClab .mpt file
            % if file not found will return empty object
            try 
                path=pwd; 
                raw_data = importdata([path '\' fileName]);
                % get table haders
                heders=raw_data.textdata(end);
                heders=heders{1};
                heders=textscan(heders,'%s','Delimiter','\t');
                temp_data=array2table(raw_data.data);
                temp_data.Properties.VariableNames=matlab.lang.makeValidName(heders{1}');
                obj.dataTable=temp_data;
                obj.E = obj.dataTable.Ewe_V;
                obj.I = obj.dataTable.x_I__mA./1000;
                obj.N_cycles = max(obj.dataTable.cycleNumber);
                speed = raw_data.textdata(35);
                speed = textscan(speed{1},'%*s %f');
                obj.speed = speed{1};
                [obj.ImaxOx,maxInd] = max(obj.I);
                obj.EmaxOx = obj.E(maxInd);
                [obj.ImaxRed,maxInd] = max(-obj.I);
                obj.EmaxRed = obj.E(maxInd);
                obj.ImaxRed = -obj.ImaxRed;
                [Ind, peakE, peakW, peakP]=findpeaks(abs(obj.I));
                peakI = obj.I(peakE); 
                peakE = obj.E(peakE);
                obj.peaks = table(peakI, peakE, peakW, peakP,Ind);
  
            catch
                obj.dataTable=[];
                obj.E = [];
                obj.I = [];
                obj.N_cycles = [];
                obj.speed = [];
                warning('You didnt entered valid File Name')
            end
                
        end
        
        function r = getExtrimom(obj)
            % this function get the exrtimom value from CV object
            % i.e max min local peaks
            [obj.ImaxOx,maxInd] = max(obj.I);
            obj.EmaxOx = obj.E(maxInd);
            [obj.ImaxRed,maxInd] = max(-obj.I);
            obj.EmaxRed = obj.E(maxInd);
            obj.ImaxRed = -obj.ImaxRed;
            [Ind, peakE, peakW, peakP]=findpeaks(abs(obj.I));
            peakI = obj.I(peakE); 
            peakE = obj.E(peakE);
            obj.peaks = table(peakI, peakE, peakW, peakP,Ind);
            r = obj;
        end
        
        function varargout=plot(obj,varargin)
            % overload plot function and enable the use of EClab_CV class
            % can get an array of objects and name value pairs and plot all
            % the objects
            pairs='';
            for i=2:nargin % use for name value pairs
                if ~isnumeric(varargin{i-1})
                    pairs=[pairs ',''' varargin{i-1} ''''];
                else
                   pairs=[pairs ',' num2str(varargin{i-1})] ;
                end
            end
             for i=1:numel(obj) % use to plot array of objects
                    if nargout==1
                        varargout{1}(i)=EClab_CV.plot_(obj(i),pairs);
                    else
                        EClab_CV.plot_(obj(i),pairs);
                    end
                    hold on
             end
             hold off
        end
         
        function r = cycleSplit(obj,varargin)
            % split the object into seprate cycles,
            % the function can get cycle index array and return only thouse
            % cycles
            % if enterd array of CVs will return object matrix in the
            % folowing fotmat:
            % CV_1_cycle_1 | CV_2_cycle_1 | ... | CV_n_cycle_1
            % CV_1_cycle_2 | CV_2_cycle_2 | ... | CV_n_cycle_2
            %       ...    |    ...       | ... |   ....
            % CV_1_cycle_m | CV_2_cycle_m | ... | CV_n_cycle_m
    
            % understand if user specified the cycles index
            if nargin==1
                cycles = 1:obj(1).N_cycles;
            else 
                cycles = varargin{1};
            end
            
            tempArrey = EClab_CV.CV_alloc(numel(cycles),numel(obj));
            % go throgh all the CV objects and split to cycles
            for j=1:numel(obj)
                for i=1:numel(cycles)
                    tempArrey(i,j) = obj(j);
                    tempArrey(i,j).N_cycles = 1;
                    tempArrey(i,j).E = obj(j).E(obj(j).dataTable.cycleNumber==cycles(i));
                    tempArrey(i,j).I = obj(j).I(obj(j).dataTable.cycleNumber==cycles(i));
                    tempArrey(i,j).dataTable = obj(j).dataTable(find(obj(j).dataTable.cycleNumber==cycles(i)),:);
                    tempArrey(i,j).dataTable.cycleNumber = ones(size(tempArrey(i,j).dataTable.cycleNumber));
                    tempArrey(i,j) = getExtrimom(tempArrey(i,j));

                end
            end
            r = tempArrey;
        end
        
        
        
        function r = minus(a,b)
        % this function used to to minuse operator on CV object
        % can use for an array of object element wise
            if numel(a)==numel(b) % if a abd b are arrays make sure that they are in the same length
                for i=1:numel(a)
                    r(i) = EClab_CV.minus_(a(i),b(i));
                end
            else
                warning('vectors must be in the same size')
            end
        end
        
        function r = getOx(a)
           % This function return the Ox part of CV as CV object
           r=a;
           N=numel(a);
           for i=1:N
               r(i).I=a(i).I(a(i).dataTable.ox_red==1);
               r(i).E=a(i).E(a(i).dataTable.ox_red==1);
               r(i)=getExtrimom(r(i));
               r(i).dataTable=a(i).dataTable(a(i).dataTable.ox_red==1,:);
           end
            
        end
        
        function r = getRed(a)
           % This function return the Ox part of CV as CV object
           r=a;
           N=numel(a);
           for i=1:N
               r(i).I=a(i).I(a(i).dataTable.ox_red==0);
               r(i).E=a(i).E(a(i).dataTable.ox_red==0);
               r(i)=getExtrimom(r(i));
               r(i).dataTable=a(i).dataTable(a(i).dataTable.ox_red==0,:);
               
           end
            
        end
        
         function scatter(obj,varargin)
        % define the scatter function for CV object
            pairs='';
            for i=2:nargin % use for name value pairs
                if ~isnumeric(varargin{i-1})
                    pairs=[pairs ',''' varargin{i-1} ''''];
                else
                   pairs=[pairs ',' num2str(varargin{i-1})] ;
                end
            end
            
            for i=1:numel(obj)
                eval(['scatter(obj(i).E,obj(i).I' pairs ')'])
                hold on
            end
            xlabel('E [V]')
            ylabel('I [A]')
            hold off
         end
        
         function array=Iarray(a)
             % This function gets CV array and convert it to I matrix
             % each coloumn is one CV
             a=evenup(a);
             N=numel(a);
             M=length(a(1).I);
             array=zeros(M,N);
             for i=1:N
                 array(:,i)=a(i).I;
             end
         end
          
         function a=evenup(b)
            % this function interpulate an arrey so the data will be in the
            % same length and points
            
            % pre alocation
            warning('off','all')
            a(size(b))=EClab_CV();
            warning('on','all')
            % get the min I samples from the array of objects
            minSamples= length(b(1).I);
            minInd=1;
            N=numel(b);
            
            for i=2:N
                if length(b(i).I)<minSamples
                    minSamples=length(b(i).I);
                    minInd=i;
                end
            end
            x=linspace(0,1,minSamples);
           
            
            for i=1:N
                
                    current_length=length(b(i).E);
                    % interpulate to even spaced vector
                    tempI=interp1(linspace(0,1,current_length),b(i).I,x);
                    tempE=interp1(linspace(0,1,current_length),b(i).E,x);
                    % update the fields
                    a(i)=b(minInd);
                    a(i).E=tempE';
                    a(i).I=tempI';
                    a(i)=getExtrimom(a(i));
                
                
            end
            
            
            
         end
        
         function r = mean(obj)
             % get obj array and return the mean current vactor
            Imat=Iarray(obj);
             r=mean(Imat,2);
         end
         
         function r = std(obj)
             % get obj array and return the std current vactor
             Imat=Iarray(obj);
             r=std(Imat,0,2);
         end
        
         function rsd = RSD(a)
             % this function gets CV coulomn arrey and return the mean RSD between
             % the signals
             % if a is a matrix it will calculate the RSD seperatly for
             % each column
             if size(a,2)>1
                 rsd=zeros(1,size(a,2));
                 for i=1:size(a,2)
                     rsd(i)=RSD(a(:,i));
                 end
             else
                 
                 rsd = mean(std(a)./abs(mean(a)))*100;
             end
         end
         
         function [base, slope,intercept]=half_cycle_baseline(single_ox)
        % This function gets half a CV cycle as an input and find the base line
            tempE=single_ox.E;
            tempI=single_ox.I;
        % this part compansate on the fact that cycle could start from the
        % middel of the Red\Ox part
            red=single_ox.dataTable.ox_red;
            if sum(red)==0
                [tempE,ind]=sort(tempE,'descend');
            else
                [tempE,ind]=sort(tempE);
            end
            tempI=tempI(ind);
            % smooth the diffrential
            Idiff=filter(ones(1,10),10,diff(tempI)./diff(tempE));
            % find 20 samples after the first local minima 
            flag=0;
            i=5;
            ind=5;
            while flag<20
                i=i+1;
                if Idiff(ind)>Idiff(i)
                    ind=i;
                    flag=0;
                else
                    flag=flag+1;
                end
            end

            IndRange=ind-4:ind+150;
          % plot(single_ox.E,single_ox.I);hold on;plot(tempE(IndRange),tempI(IndRange))
            % fit a slope and intercept
            fitObj=fitlm(tempE(IndRange),tempI(IndRange));
            slope=fitObj.Coefficients.Estimate(2);
            intercept=fitObj.Coefficients.Estimate(1);
            base=intercept+single_ox.E.*slope;
            
            
        end
         
         function [baseMat,varargout]=Baseline_Est(obj)
             % this function gets CV and return CV matrix in this shape
             % ||cycle 1|cycle 1 ox baseline | cycle 1 Red baseline ||
             % ||cycle 2|cycle 2 ox baseline | cycle 2 Red baseline ||
             % ||cycle 3|cycle 3 ox baseline | cycle 3 Red baseline ||
             % [baseMat,NoBase]= Baseline_Est(obj) will give a signal
             % without the baseline

           cycled=cycleSplit(obj);
           for i=1:length(cycled)
               oxpart(i)=getOx(cycled(i));
               [oxpart(i).I,slopeOx,inteceptOx]=half_cycle_baseline(oxpart(i));
               redpart(i)=getRed(cycled(i));
               [redpart(i).I,slopeRed,inteceptRed]=half_cycle_baseline(redpart(i));
               NoBase(i)=cycled(i);
               Base(i)=cycled(i);
               NoBase(i).I=cycled(i).I-(inteceptOx+slopeOx.*cycled(i).E).*(cycled(i).dataTable.ox_red==1)-(inteceptRed+slopeRed.*cycled(i).E).*(cycled(i).dataTable.ox_red==0);
               NoBase(i)=getExtrimom(NoBase(i));
               Base(i).I=(inteceptOx+slopeOx.*cycled(i).E).*(cycled(i).dataTable.ox_red==1)+(inteceptRed+slopeRed.*cycled(i).E).*(cycled(i).dataTable.ox_red==0);
            
           end
           baseMat=[cycled,oxpart',redpart'];
           switch nargout
                case 2
                    varargout{1}=NoBase;
                case 3
                    varargout{2}=Base;
            end
         end

        function [E,I]=getPeaks(CV,N)
            % this function returns peak E and I of N most prominant peaks
            sorted_peaks=sort(CV.peaks.peakP,'descend');
            % get the N biggest peaks
            p=sorted_peaks(N);
            if numel(CV)==1 % if for a singal object
                % get the N bigest peaks wich is not in the end or
                % baggining of the signal
                E=CV.peaks.peakE(CV.peaks.peakP>=p & CV.peaks.peakE<(max(CV.peaks.peakE)-0.002) & CV.peaks.peakE>(min(CV.peaks.peakE)+0.002));
                I=CV.peaks.peakI(CV.peaks.peakP>=p & CV.peaks.peakE<(max(CV.peaks.peakE)-0.002) & CV.peaks.peakE>(min(CV.peaks.peakE)+0.002));
                E(E==0)=NaN;
                I(I==0)=NaN;
            else % if for array of CVs
                E=zeros(50,numel(CV));
                I=zeros(50,numel(CV));
                for i=1:numel(CV)
                    [tempE,tempI]=getPeaks(CV(i),N);
                    E(1:length(tempE),i)=tempE;
                    I(1:length(tempE),i)=tempI;
                end
                E( all( ~any(E,2), 2 ), : ) = [];
                I( all( ~any(I,2),2), : ) = [];
            end
        end
        
        function [Base,NoBase]=BaseSubtract(obj)
            if numel(obj)>1
                for i=1:numel(obj)
                      [Base(i),NoBase(i)]=BaseSubtract(obj(i))
                end
            else
            cycled=cycleSplit(obj);
                for i=1:length(cycled)
                   oxpart=getOx(cycled(i));
                   [~,slopeOx,inteceptOx]=half_cycle_baseline(oxpart);
                   redpart=getRed(cycled(i));
                   [~,slopeRed,inteceptRed]=half_cycle_baseline(redpart);
                   cycleNoBase(i)=cycled(i);
                   cycleBase(i)=cycled(i);
                   cycleNoBase(i).I=cycled(i).I-(inteceptOx+slopeOx.*cycled(i).E).*(cycled(i).dataTable.ox_red==1)-(inteceptRed+slopeRed.*cycled(i).E).*(cycled(i).dataTable.ox_red==0);
                   cycleBase(i).I=(inteceptOx+slopeOx.*cycled(i).E).*(cycled(i).dataTable.ox_red==1)+(inteceptRed+slopeRed.*cycled(i).E).*(cycled(i).dataTable.ox_red==0);
               end
               Base=obj;
               Base.I=[];
               NoBase=obj;
               NoBase.I=[];
               for i=1:numel(cycleBase)
                    Base.I=[Base.I;cycleBase(i).I];
                    NoBase.I=[NoBase.I;cycleNoBase(i).I];

               end
               Base=getExtrimom(Base);
               NoBase=getExtrimom(NoBase);

            end
        
           
        end
  
         
    end
end