%
%
% This example script performs the analysis included in Supplementary Note 6.
% We simulate guided wave propagation due to a line source incident on a 
% spherical isotropic layer bounded above by air and below by water. This 
% provides a model for guided wave propagation in the curved geometry
% of the cornea. This script shows that along the line normal to the center
% of the line source, the XT and FK plots for the curved surface are nearly
% identical to those of a flat surface. This suggests that a flat plate model
% of the cornea captures the essential behavior and curvature can be neglected.
%
% The outputs for this script are a movie of wave propagation over a curved 
% surface (Supplementary Video 1), and XT and FK plots mimicking an OCE 
% measurement (Supplemental Figure S6.2). 
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%


%
% Load OnScale results for a point source on a curved surface
% For the sake of this example, we've provided the OnScale results
% already loaded into a *.mat file with the correct data structure.
% If you run the OnScale model point_source_curved.flxinp, you can
% generate this data structure by calling:
%
% data = load_onscale_flexdata(folder_containing_results, 'psc_results', 1);
%
% Here, folder_containing_results is the path to your OnScale results,
% each of which is a numbered *.flxdato file with the prefix 'psc_results'.
%
load('example_onscale_results_point_source_curved.mat');
dx = data.x(2) - data.x(1);

%
% Generate grid for domain
% sph_radius is the radius used in creating the curved domain in OnScale.
% Here, it approximates the corneal radius of curvature.
% npts is the number of nodes in the surface mesh.
% Points are distributed using a golden spiral method
%
sph_radius = 0.0065;
npts = 40000;
idx = ((0:(npts-1)) + 0.5)';
phi = acos(1 - 2*idx/npts);
theta = pi * (1 + sqrt(5))*idx;
smesh.x = sph_radius*cos(theta).*sin(phi);
smesh.y = sph_radius*sin(theta).*sin(phi);
smesh.z = sph_radius*cos(phi);

% 
% Due to symmetry, it is sufficient to consider only 1/4 of the hemisphere. 
%
r = sqrt(smesh.x.^2 + smesh.z.^2);
mask = (r <= 0.006) & (smesh.y >= 0) & (smesh.x >= 0);
smesh.x = smesh.x(mask);
smesh.y = sph_radius - smesh.y(mask);
smesh.z = smesh.z(mask);

%
% Compute a Delaunay triangulation of the mesh - this is not needed for
% computing solutions, but is useful for plotting results.
%
tri = delaunayTriangulation(smesh.x, smesh.y, smesh.z);
smesh.elmcon = tri.convexHull;

%
% Generate source points for the line source
% lensrc is the length of the line source, here 6 mm to approximate our
% acoustic micro-tapping (AuT) transducer's acoustic field. The number
% of source points is set so that it uses roughly the same spatial sampling
% as the OnScale model. The line source points are mapped over the spherical
% surface.
%
lensrc = 0.006;
nsrc = round(lensrc/dx) + 1;
thsrc = linspace(-asin(lensrc/2/sph_radius), asin(lensrc/2/sph_radius), nsrc)';
xsrc = zeros(nsrc, 1);
ysrc = sph_radius - sph_radius*cos(thsrc);
zsrc = sph_radius*sin(thsrc);


%
% Show the mesh
%
figure
set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
tplt = trisurf(smesh.elmcon, smesh.x*1000, smesh.z*1000, smesh.y*1000);
tplt.FaceColor = [0.8, 0.8, 0.8];
axis equal
hold on
plot3(xsrc*1000, zsrc*1000, ysrc*1000, '.r')
set(gca, 'ZDir', 'reverse')
set(gca, 'XDir', 'reverse')
plot3([5, 4], [5, 5], [1, 1]*max(smesh.y)*1000, '-k', 'LineWidth', 1)
plot3([5, 5], [5, 4], [1, 1]*max(smesh.y)*1000, '-k', 'LineWidth', 1)
plot3([5, 5], [5, 5], [1, 1]*max(smesh.y)*1000 - [0, 1], '-k', 'LineWidth', 1)
axis off
print('half_globe_domain_preview_fine', '-dpng', '-r300')


%
% Solve the convolution of the point source solution with the line source distribution
% 
dataout = onscale_point_source_curved_integral(data, xsrc, ysrc, zsrc, smesh, sph_radius);


% 
% Create a movie showing the wave propagation results, similar to Supplementary Video 1
%
filename = 'propagation_movie';
vid = VideoWriter(filename, 'MPEG-4');
vid.FrameRate = 20;
vid.Quality = 90;

hfig = figure;
set(gcf, 'Color', [0, 0, 0], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
axes
set(gca, 'Color', [0, 0, 0]);
vnorm = max(max(abs(dataout.v(data.x < 0.01, dataout.t < 0.002))));
open(vid);
for j = 1:length(dataout.t)
    if dataout.t(j) > 0.002
        break
    end
    cla
    pp = patch('Faces', dataout.smesh.elmcon, 'Vertices', [dataout.smesh.x, dataout.smesh.z, dataout.smesh.y]*1000, 'FaceVertexCData', dataout.v(:,j)/vnorm, 'FaceColor', 'interp', 'LineStyle', 'none');
    axis equal
    view(3)
    set(gca, 'XDir', 'reverse')
    set(gca, 'ZDir', 'reverse')
    xticks([])
    yticks([])
    zticks([])
    tt = title(['Time = ', num2str(dataout.t(j)*1000, '%1.3f'), ' ms']);
    tt.FontWeight = 'normal';
    tt.Color = 'w';
    caxis([-0.5, 0.5])
    cb = colorbar;
    ylabel(cb, 'Surface Vertical Velocity (a.u.)')
    hold on
    plot3([5, 4], [5, 5], [1, 1]*max(smesh.y)*1000, '-w', 'LineWidth', 1)
    plot3([5, 5], [5, 4], [1, 1]*max(smesh.y)*1000, '-w', 'LineWidth', 1)
    plot3([5, 5], [5, 5], [1, 1]*max(smesh.y)*1000 - [0, 1], '-w', 'LineWidth', 1)
    xlim([0, 6.5])
    ylim([-6.5, 6.5])
    zlim([-0.5, max(smesh.y)*1000 + 0.5])
    cb.Label.Color = [1,1,1];
    cb.Label.FontSize = 10;
    set(gca, 'XColor', 'k')
    set(gca, 'YColor', 'k')
    set(gca, 'ZColor', 'k')
    cb.Color = 'w';
    t2 = text(6.3, 5.5, 3.1, '1 mm');
    t2.Color = 'w';
    
    frame = getframe(hfig);
    writeVideo(vid, frame);

end
close(vid)


%
% In the following section, we re-compute the solution only along the OCT imaging plane
% This is a line normal to the line source center. This mirrors our OCE measurement, and so
% we want to compare XT plots and FK spectra acquired along this line. This mimics the plots
% in Supplementary Figure S6.2.
%

%
% Generate grid for domain
%
sph_radius = 0.0065;
smesh.x = data.xs(:);
smesh.y = data.ys(:);
smesh.z = zeros(size(smesh.x));

%
% Generate source points
%
lensrc = 0.006;
nsrc = round(lensrc/dx) + 1;
thsrc = linspace(-asin(lensrc/2/sph_radius), asin(lensrc/2/sph_radius), nsrc)';
xsrc = zeros(nsrc, 1);
ysrc = sph_radius - sph_radius*cos(thsrc);
zsrc = sph_radius*sin(thsrc);

%
% Solve
%
dataout2 = onscale_point_source_curved_integral(data, xsrc, ysrc, zsrc, smesh, sph_radius);
vnorm2 = max(max(abs(dataout2.v(data.x < 0.01, dataout2.t < 0.002))));

%
% Display XT plot for the curved solution
%
ds = sqrt(diff(data.xs).^2 + diff(data.ys).^2);
ss = [0, cumsum(ds)];
figure
pcolor(data.t*1000, ss*1000, dataout2.v/vnorm2)
shading flat
xlabel('Time (ms)')
cb = colorbar;
ylabel(cb, 'Surface Vertical Velocity (a.u.)');
caxis([-0.5, 0.5])
xlim([0, 2])
ylim([0, 10])
ylabel('Lateral Position (mm)')


%
% Compute power spectrum for FK plot
%
dt = data.t(2) - data.t(1);
nfft = 4096;
[~, idx] = min(abs(data.x - 0.01));
xrange = 1:idx;
[~, idx] = min(abs(data.t - 0.002));
trange = 1:idx;
[VdB, f, k] = xttools_power_spectrum(dataout2.v(xrange, trange), dt, dx, nfft);


%
% Display FK plot
%
figure; 
pcolor(f/1000, k/1000, VdB); 
shading flat; 
colormap(hot); 
caxis([-20, 0])
cb = colorbar;
ylabel(cb, 'Power Spectrum (dB)')
xlim([0, 8])
ylim([0, 2])
xlabel('Frequency (kHz)')
ylabel('Wavenumber (1/mm)')

