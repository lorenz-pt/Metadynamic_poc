%% Plot of the charge value
% The program print a figure with two subplot, one showing the evolution of
% the charge (the one on the left) and the other showing the histogram of
% occurrences.

%Use this for loading from file;
%Q = load('.txt');
%if you have the data saved on the workspace;
Q = Q_min;
%Sweeps during the simulation;
sweeps = linspace(1,length(Q),length(Q));
%set the title;
sgtitle('$\beta\chi\hbar = 0.5$, sweeps = 1e5, nlattice = 300','interpreter','latex','Fontsize',16);
%Two subplot for showing the evolution of the charge with an histogram and
%with the evolution during the simulation;
subplot('Position',[.1,0.1,.6,.8]);

h1 = plot(sweeps,Q);
h1.Color = [.8,.1,.2];

%Label for the plot of the charge;
ylabel('Q','interpreter','latex','Fontsize',16);
xlabel('Sweeps','interpreter','latex','Fontsize',16);

subplot('Position',[.75,.1,.2,.8]);
%histogram of the charge;
h2 = histogram(Q);

%Colors of the histogram;
h2.FaceColor = [.8,.1,.2];
h2.EdgeColor = [1,0,0];
%set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
box off;
%Rotated view;
set(gca, 'view',[90 -90]);
set(gcf, 'PaperSize',[20,16]);

