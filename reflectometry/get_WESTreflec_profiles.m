function [ imas_edge_density,imas_core_density_REDini, imas_core_density_INTini,tr,ic,ia,profil] ...
    = get_WESTreflec_profiles(shot,profil)
% [ imas_edge_ne,imas_core_density_REDinit, imas_core_density_REDinit] = get_WESTreflec_profile(shot )
% get the density profiles from WEST edge and core reflectometers
% imas_edge_ne: edge density profile from DREFRAP
% imas_core_density_REDini: core density profiles from DREFLUC initilized on DREFRAP profiles 
% imas_core_density_INTini: core density profiles initilized on interferometry (IMAS core_profiles)
% the fonction recovers the edge profiles used for the initilization 
% the complete profiles (position (r,psi) and density n_e
% are stored as additionnal fields on the root of the IMAS structure (see future IMAS version) 
% 

% R Sabot (06/2021)
% verification performed on pulse 55987 (07/06/21) 

if nargin>1 && isstruct(profil)
    noprofil_in=false;
else
    noprofil_in=true;
end
[tr,ic,ia]=deal([]);
Nout=128;
if length(shot)==1 && shot>50000
    try 
        imas_edge_density=imas_west_get(shot,'reflectometer_profile',0,0);
    catch
        imas_edge_density.time=[];
    end
    try 
        imas_core_density_REDini=imas_west_get(shot,'reflectometer_profile',0,1);
    catch
        imas_core_density_REDini.time=[];
    end
    try
        imas_core_density_INTini=imas_west_get(shot,'reflectometer_profile',0,2);
    catch
        imas_core_density_INTini.time=[];
    end
    
    if isempty(imas_core_density_INTini.time) 
        fprintf('\nNo DREFLUC data for %d ',shot )
        imas_core_density_INTini.n_e.data=[];
        imas_core_density_INTini.position.psi.data=[];
        imas_core_density_INTini.position.r.data=[];
    end
    if isempty(imas_core_density_REDini.time) 
        disp('\nNo imas data with edge reflectometry initialisation')
        imas_core_density_REDini.n_e.data=[];
        imas_core_density_REDini.position.psi.data=[];
        imas_core_density_REDini.position.r.data=[];
    end
%     if isempty(imas_core_density_REDini.time) && isempty(imas_core_density_INTini.time) 
%         return;
%     end
    N_INT=length(imas_core_density_INTini.time);
    N_RED=length(imas_core_density_REDini.time);
    N_FRAP=length(imas_edge_density.time);
        
    if shot<56930
        if noprofil_in 
            profil=imas_west_get(shot, 'core_profiles',0, 0,'imas_public','west');
        end
        [tr,ia,ic]=unique([imas_core_density_INTini.time;imas_core_density_REDini.time;...
            imas_edge_density.time]);
        ic_INT=ic(1:N_INT);
        ic_RED=ic(N_INT+1:N_INT+N_RED);
        ic_FRAP=ic(N_INT+N_RED+1:end);
        equi=equi_freq_local(shot,tr,ic_INT,profil);
        Rcore_INT=imas_core_density_INTini.channel{1}.position.r.data;
        necore_INT=imas_core_density_INTini.channel{1}.n_e.data;
        Rcore_RED=imas_core_density_REDini.channel{1}.position.r.data;
        necore_RED=imas_core_density_REDini.channel{1}.n_e.data;
%         [Redge_INT,needge_INT]=deal(zeros(Nedge,N_INT));
        imas_core_density_INTini.position.r.data=zeros(Nout,N_INT);
        imas_core_density_INTini.n_e.data=zeros(Nout,N_INT);
%         [Redge_RED,needge_RED]=deal(zeros(Nedge,N_RED));
        
        for ii=1:N_INT
            iiequi=ic(ii);
%            Rtemp=linspace(Rcore_INT(end,iiequi),equi.R(end),Nedge+1)';
%            Redge_INT(:,ii)=Rtemp(2:end);

            imas_core_density_INTini.position.r.data(:,ii)=linspace(Rcore_INT(1,ii),equi.R(end),Nout)';
            indR=equi.R>Rcore_INT(end,ii);
%            needge_INT(:,ii)=interp1(equi.R,equi.ne(iiequi,:),Redge_INT(:,ii));
            imas_core_density_INTini.n_e.data(:,ii)=interp1([Rcore_INT(:,ii);equi.R(indR)'],...
                [necore_INT(:,ii);equi.ne(iiequi,indR)'],imas_core_density_INTini.position.r.data(:,ii),'linear');
%             imas_core_density_INTini.n_e.data(:,ii)=interp1([Rcore_INT(:,ii);Redge_INT(:,ii)],...
%                 [necore_INT(:,ii);needge_INT(:,ii)],imas_core_density_INTini.position.r.data(:,ii),'linear');
%             imas_core_density_INTini.position.psi.data(:,ii)=interp1(equi.R,equi.psi(iiequi,:),...
%                 imas_core_density_INTini.position.r.data(:,ii),'linear');                      
            imas_core_density_INTini.position.rho_pol_norm.data(:,ii)=interp1(equi.R,equi.rho_pol_norm(iiequi,:),...
                imas_core_density_INTini.position.r.data(:,ii),'linear');                      
            imas_core_density_INTini.position.rho_tor_norm.data(:,ii)=interp1(equi.R,equi.rho_tor_norm(iiequi,:),...
                imas_core_density_INTini.position.r.data(:,ii),'linear');                      
            imas_core_density_INTini.position.theta.data(:,ii)=interp1(equi.R,equi.theta(iiequi,:),...
                imas_core_density_INTini.position.r.data(:,ii),'linear');                      
        end
        prev_frap=1;
        for ii=1:N_RED
            iiequi=ic_RED(ii);
%             Rtemp=linspace(Rcore_RED(end,ii),equi.R(end),Nedge+1)';
%             Redge_RED(:,ii)=Rtemp(2:end);
            imas_core_density_REDini.position.r.data(:,ii)=linspace(Rcore_RED(1,ii),equi.R(end),Nout)';
            
            [~,indice]=min(abs(imas_edge_density.time(prev_frap:end)-imas_core_density_REDini.time(ii)));
            ind_frap=indice+prev_frap-1;

            [RXs,NEXs]=reordonne(imas_edge_density.channel{1}.position.r.data(:,ind_frap),...
                imas_edge_density.channel{1}.n_e.data(:,ind_frap));
            indRX1=RXs>Rcore_RED(end,ii);
            imas_core_density_REDini.n_e.data(:,ii)=interp1([Rcore_RED(:,ii);RXs(indRX1)],...
                [necore_RED(:,ii);NEXs(indRX1)],...
                imas_core_density_REDini.position.r.data(:,ii),'linear',0);
            
%             imas_core_density_REDini.n_e.data(:,ii)=interp1([Rcore_RED(:,ii);Redge_RED(:,ii)],...
%                 [necore_RED(:,ii);needge_RED(:,ii)],imas_core_density_REDini.position.r.data(:,ii),'linear');
            imas_core_density_REDini.position.rho_pol_norm.data(:,ii)=interp1(equi.R,equi.rho_pol_norm(iiequi,:),...
                imas_core_density_REDini.position.r.data(:,ii),'linear');
            imas_core_density_REDini.position.rho_tor_norm.data(:,ii)=interp1(equi.R,equi.rho_tor_norm(iiequi,:),...
                imas_core_density_REDini.position.r.data(:,ii),'linear');
            imas_core_density_REDini.position.theta.data(:,ii)=interp1(equi.R,equi.theta(iiequi,:),...
                imas_core_density_REDini.position.r.data(:,ii),'linear');
            

        end
 % to complete DREFRAP
        if N_FRAP
            if ~isfield(imas_edge_density,'n_e')
                imas_edge_density.n_e.data=imas_edge_density.channel{1}.n_e.data;
            end
            if ~isfield(imas_edge_density,'position')
                imas_edge_density.position.r.data=imas_edge_density.channel{1}.position.r.data;
            end
            imas_edge_density.position.rho_tor_norm.data=NaN(size(imas_edge_density.position.r.data));
            imas_edge_density.position.rho_pol_norm.data=NaN(size(imas_edge_density.position.r.data));
            for ii=1:N_FRAP
                iiequi=ic_FRAP(ii);
                imas_edge_density.position.rho_tor_norm.data(:,ii)=interp1(equi.R,equi.rho_tor_norm(iiequi,:),...
                    imas_edge_density.position.r.data(:,ii),'linear');
                imas_edge_density.position.rho_pol_norm.data(:,ii)=interp1(equi.R,equi.rho_pol_norm(iiequi,:),...
                    imas_edge_density.position.r.data(:,ii),'linear');
            end
        else
            imas_edge_density.n_e.data=[];
            imas_edge_density.position.r.data=[];
            imas_edge_density.position.rho_pol_norm.data=[];
            imas_edge_density.position.rho_tor_norm.data=[];
        end
    end
    
else
    warning('shot must be a WEST pulse, ie an integer > 50000');
%     listeshot=shot;
%     for ii=1:numel(listeshot)
%         shot=listeshot(ii);
%         imas_edge_density{ii}=imas_west_get(shot,'reflectometer_profile',0,0);
%         imas_core_density_REDini{ii}=imas_west_get(shot,'reflectometer_profile',0,1);
%         imas_core_density_INTini{ii}=imas_west_get(shot,'reflectometer_profile',0,2);
%     end
end

end

function equi=equi_freq_local(pulse,tr,ic_INT,profil)
% equi=equi_freq(pulse,tr,profil)

% tr : time reflectometry / time must be absolute (w/o ignitron)
% travail sur equi.fpe
% Hyp: one reflectometer or all reflectometers at (Z=0, phi=0)
% R. Sabot
% 02/07/2020: correction sur rho outside rho=1:
% sqrt only on psi, not on prolongation outside LCFS
% 
%global ReLFS ReHFS Rstep lambda_nerho lambda_Terho
% global data
ReLFS=3.1; 
ReHFS=2.0; % ???
Rstep=0.005;
lambda_nerho=0.04; %ie ~ 2 cm


interpol_method='linear';

equi.R=ReHFS:Rstep:ReLFS;
equi.lR=length(equi.R);
equi.Z=zeros(size(equi.R));
equi.phi=zeros(size(equi.R));

equi.time=tr;
equi.lt=length(tr);
occur_equi_Ini=0;
equi.occur_comment='NICE only';% & polarimetry';
comment_occ={'NICE only','NICE & polarimetry'};
try
    out = imas_west_get(pulse, 'equilibrium',0,occur_equi_Ini);
    if isempty(out) || isempty(out.time)
       occur_equi=mod(occur_equi_Ini+1,2); 
    else
        occur_equi=occur_equi_Ini;
    end
catch
    occur_equi=mod(occur_equi_Ini+1,2);
end
if occur_equi~=occur_equi_Ini;
    out = imas_west_get(pulse, 'equilibrium',0,occur_equi);
    if isempty(out) || isempty(out.time)
        error('TRCO:equi_freq','No equilibrium');
    end
end    
equi.occur_equi=occur_equi;
equi.occur_comment=comment_occ{1+occur_equi};%'NICE & polarimetry';

equi.equilibrium=out;

% equi.Btotg=equimap_get(pulse, equi.time, equi.R, equi.phi, equi.Z, 'b_field_norm',false,0,occur_equi);
equi.theta=equimap_get(pulse, equi.time, equi.R, equi.phi, equi.Z, 'theta',false,0,occur_equi);
equi.psi=equimap_get(pulse, equi.time, equi.R, equi.phi, equi.Z, 'psi',false,0,occur_equi);
equi.rho_pol_norm=equimap_get(pulse, equi.time, equi.R, equi.phi, equi.Z, 'rho_pol_norm',false,0,occur_equi);
equi.rho_tor_norm=equimap_get(pulse, equi.time, equi.R, equi.phi, equi.Z, 'rho_tor_norm',false,0,occur_equi);
%equi.rho_pol_normS=equi.rho_pol_norm.*sign(cos(equi.theta));
equi.indR0=min(equi.rho_pol_norm,[],1);

equi.ne=zeros(equi.lt,equi.lR);

% interpolation at T_INT time to get the edge profile

dtms=round(median(diff(profil.time))*1000);
dt=0.001*dtms;
profil.dtms=dtms;



tr_INT=tr(ic_INT);
N_INT=length(ic_INT);
%N_INT=length(Ic_INT);
% if strcmp(interpol_method,'linear')
% linear interpolation
% we can select the profil.times that are around reflectometry time
ind_trinf=1+floor((tr_INT-profil.time(1))/dt);
ind_trinf(ind_trinf>=length(profil.time))=length(profil.time)-1;
indpb=find(tr_INT(:)-profil.time(ind_trinf)<=0);
if ~isempty(indpb)
    ind_trinf=[ind_trinf;max(1,ind_trinf(indpb)-1)];
end
ind_trinfU=unique(ind_trinf);
I_profS_INT=unique([ind_trinfU;ind_trinfU+1]);
% else
%     I_profS_INT=(1+floor((tr_INT(1)-profil.time(1))/dt)):(length(profil.time)-floor((profil.time(end)-tr_INT(end))/dt));
% end
% I_profS=ic_INT(I_profS_INT);

 Nof_tprof=length(I_profS_INT);

%rhomax=; % to get rho_pol grid up to 1.15
lpsi=length(profil.profiles_1d{1}.grid.psi);

% if lpsi==101
%     psing=0:0.01:psinmax;
% else
%    psing=[linspace(0,1,lpsi),1+1/(lpsi-1):1/(lpsi-1):psinmax];
% end
drho=0.01;% the profile grid is a grid on psi 0:0.01:1;
% In 10 ms, the position of the flux surfaces do not evolve much 
% the psi grid is supposed to be fixed (in R) 
% and one can perform a simple the time interpolation on this grid

equi.psi=linspace(0,1,lpsi);
rho_pol_ng=[sqrt(equi.psi),1+drho:drho:1+2.5*lambda_nerho,1+3*lambda_nerho:lambda_nerho/2:1+5*lambda_nerho,...
    1+6*lambda_nerho:lambda_nerho:1+8*lambda_nerho];

lrhog=length(rho_pol_ng);


ne_psi=zeros(lpsi,Nof_tprof);
% Te_psi=zeros(lpsi,Nof_tprof);

for ii=1:Nof_tprof
    indprof=I_profS_INT(ii);
    ne_psi(:,ii)=profil.profiles_1d{indprof}.electrons.density;
%     Te_psi(:,ii)=profil.profiles_1d{indprof}.electrons.temperature;
end
% replace Nan Te-0 => no relativistic correction
% pb_Te=isnan(Te_psi(:));
% if any(pb_Te)
%     Te_psi(pb_Te)=0;
% end
% if unexpected ne = NaN => ne =0
pb_ne=isnan(ne_psi(:));
if any(pb_ne)
    ne_psi(pb_ne)=0;
end

ne_psi_tr_INT=interp1(profil.time(I_profS_INT),transpose(ne_psi),tr_INT,interpol_method,0);
% Te_psi_tr=interp1(profil.time(I_profS),transpose(Te_psi),tr,interpol_method,0);

% remove negative value

% negat=find(ne_psi_tr<0);
% if ~isempty(negat)
%     ne_psi_tr(negat)=0;
% end
% neg_Te=(Te_psi_tr(:)<0 );
% if any(neg_Te)
%     Te_psi_tr(neg_Te)=0;
% end
neg_ne=(ne_psi_tr_INT(:)<0 | isnan(ne_psi_tr_INT(:)));
if any(neg_ne) 
    ne_psi_tr_INT(neg_ne)=0;
end
% if any(isnan(Te_psi_tr(:)))
%     keyboard;
% end

ne_psi_tr_INT(:,lpsi+1:lrhog)=bsxfun(@times,ne_psi_tr_INT(:,lpsi),exp((1-rho_pol_ng(lpsi+1:lrhog))/lambda_nerho));
%Te_psi_tr(:,lpsi+1:lrhog)=bsxfun(@times,Te_psi_tr(:,lpsi),exp((1-rho_pol_ng(lpsi+1:lrhog))/lambda_Terho));

% equi.Te=zeros(equi.lt,equi.lR);
% equi.fce=zeros(equi.lt,equi.lR);
% equi.fcxh=zeros(equi.lt,equi.lR);
% equi.fpe=zeros(equi.lt,equi.lR);
% equi.cor_relat=zeros(equi.lt,equi.lR);

% since global interpolation 
%equi.ne=interp1(rho_pol_ng,transpose(ne_psi_tr),transpose(equi.rho_pol_norm));
% doesn't work fast, better to perform the interpolation in the loop
% interp1(rho_pol_ng,ne_psi_tr(1000,:),equi.rho_pol_norm(1000,:))
for ij=1:N_INT%equi.lt
    ii=ic_INT(ij);
    equi.ne(ii,:)=interp1(rho_pol_ng,ne_psi_tr_INT(ij,:),equi.rho_pol_norm(ii,:),interpol_method,0);
%     equi.Te(ii,:)=interp1(rho_pol_ng,Te_psi_tr(ii,:),equi.rho_pol_norm(ii,:),interpol_method,0);
end

% equi.cor_relat = sqrt((1+5e-3/511*equi.Te));
% equi.fce=27.99e9*equi.Btotg./equi.cor_relat ;
% equi.fpe=sqrt(80.61*equi.ne./equi.cor_relat);
% equi.fcxh=(equi.fce+sqrt(equi.fce.^2 + 4*equi.fpe.^2))/2;
% equi.rho2fce=1;
% equi.indrho2fce=iroundrs(equi.rho_pol_normS',equi.rho2fce);
% equi.twoFceGHz=2e-9*matINmat(equi.fce',equi.indrho2fce);
% equi.max2fce=max(equi.twoFceGHz);

% manque reflectometry 
end 

function [Rsort,ysort,Isort,Isort2]=reordonne(RI,yI,sens)
% [Rsort,ysort]=reordonne(RI,yI,sens)
% reordonne les profils de dreflec pour tjrs avoir des profils croissant en R
% fait une mayenne si 2 valeurs identiques en R
% 
% Copy from core_reflectometer data processing 
% R. Sabot 30/06/20: matrices OK

if min(size(RI))==1
    [Rsort,Isort]=sort(RI); 
    if nargin >1
        ysort=yI(Isort);
        egaux=find(diff(Rsort)==0);
        % il peut arriver que 2 rayons soient identiques
        % on fait alors la moyenne de ne a ces 2 rayons
        while ~isempty(egaux)
            Rsort(egaux)=Rsort(egaux+1)*(1-2*eps);
            ysort(egaux)=(ysort(egaux)+ysort(egaux+1))/2;
            ysort(egaux+1)=ysort(egaux);
            egaux=find(diff(Rsort)==0);
        end
    else
        ysort=[];
    end
    if nargout>3
        Isort2=Isort;
    end
else
    if nargin<3
        sens=1;
    elseif sens>1
        sens=2;
    else
        sens=1;
    end
    if sens==2
        RI=transpose(RI);
        yI=transpose(yI);
    end
    [nrow,ncol]=size(RI);    
    [Rsort,Isort]=sort(RI,1);
    Isort2=bsxfun(@plus,Isort,0:nrow:(ncol-1)*nrow);
    ysort=yI(Isort2);
    I0=(diff(Rsort,[],1)==0);
    while any(I0)
        tI0=[I0;false(1,ncol)];
        Rsort(tI0)=Rsort(tI0)*(1-2*eps);
        nextI0=[false(1,ncol);I0];
        vmoy=(ysort(tI0)+ysort(nextI0))/2;
        ysort(tI0)=vmoy;
        ysort(nextI0)=vmoy;
        I0=(diff(Rsort,[],1)==0);

    end
    if sens==2
       ysort=transpose(ysort);
       if nargin >2
           Isort=transpose(Isort);
       end
    end
end
        
return
end



