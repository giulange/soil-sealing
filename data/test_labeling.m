%% set dim:

% fixed:
tiledimX    = 30;
tiledimY    = 30;

% variable:
ntilesX     = 1;
ntilesY     = 8;
threshold   = 0.7;

% consequence:
NC          = (ntilesX-1)*(tiledimX-1) + tiledimX -2;
NR          = (ntilesY-1)*(tiledimY-1) + tiledimY -2;

%% save binary image

A           = rand(NR,NC);

A(A>=threshold)   = 1;
A(A<threshold)    = 0;

fid         = fopen('/home/giuliano/git/soil-sealing/data/ALL.txt','w');
for irow = 1:NR
    for icol = 1:NC
        fprintf(fid, '%d ', A(irow,icol));
    end
    fprintf(fid,'\n');
end
fclose(fid);

%% compute labeled image

% compute labeled image from within MatLab;
tic,
L           = bwlabel(A',8);
T(1) = toc;
L = L';
%% load C-code labeled image
tic
eval( ['!/home/giuliano/git/soil-sealing/Debug/soil-sealing ',num2str(tiledimX),' ',num2str(tiledimY),' ',num2str(NC),' ',num2str(NR)] );
T(2) = toc;
fid         = fopen('/home/giuliano/git/soil-sealing/data/Ccode.txt','r');
Ccode       = fscanf(fid, '%d ', [NC,NR]);
fclose(fid);

Ccode       = Ccode';

%% compare

DIFFER = Ccode - L;
idx=find(DIFFER);
UL = unique([L(idx),Ccode(idx)],'rows');
UC = unique([Ccode(idx),L(idx)],'rows');
% sum(diff(UL(:,1))==0)
% sum(diff(UC(:,1))==0)

fprintf( '%s\n', repmat('*',1,55) );
fprintf( '   The Connected Component Labeling (CCL) procedure\n' );
fprintf( '\t   # CUDA vs MatLab # \n' );
fprintf( '   Image[%d,%d], tiledim[%d,%d], ntiles[%d,%d] \n',NR,NC,tiledimX,tiledimY,ntilesX,ntilesY );
fprintf( '   differs in <<%d>> cells [norm=%.2f]\n', length(find(DIFFER)),norm(DIFFER) );
fprintf( '     > nID: %d vs %d\n',max(Ccode(:)),max(L(:)) );
fprintf( '     > %d MatLab labels have %d C-code labels\n', max(bwlabel(diff(UL(:,1)),4))-1, length(unique(UL( find(bwlabel(diff(UL(:,1)),4), 2)))) );
fprintf( '     > %d C-code labels have %d MatLab labels\n', max(bwlabel(diff(UC(:,1)),4))-1, length(unique(UC( find(bwlabel(diff(UC(:,1)),4),2)))) );
fprintf( '%s\n', repmat('*',1,55) );

%% plot

figure(1),subplot(131),gpcolor(DIFFER),colorbar,title('(C-code - MatLab) CCL differences')
figure(1),subplot(132),gpcolor(L),colorbar,title('MatLab CCL')
figure(1),subplot(133),gpcolor(Ccode),colorbar,title('C-code CCL')