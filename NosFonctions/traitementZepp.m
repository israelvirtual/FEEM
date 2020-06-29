function [eem_cor, correct, eem_filter] = cleanscan2(eem, correct, parms, ex, em, baseopt)
%syntax : [eem_cor] = cleanscan2(eem, tol, coeff, ex,
%em, baseopt)
%
% Corrects an excitation-emission matrix fluorometer scan (EEM) to remove
% remove specific optical aberrations(i.e. 1st-order Rayleigh scatter, 2nd
% order Rayleigh scatter (i.e. overtones) and their associated Raman
% scatter peaks).
% If tolerances < 1000 are specified in the 'correct' matrix, fluorescence
% values in scatter regions (i.e.. peak emission +/- tolerance in nm for
% each excitation wavelength) are excised and replaced by interpolation of
% the surrounding data using a three-dimensional Delaunay triangulation
% method. If tolerances >= 1000 are specified, values in the corresponding
% portion of the EEM will be replaced with a constant value specified by
% 'baseopt' (i.e. truncated).
% Tolerances values of 0 can be used to selectively disable correction of
% the scatter peaks
%
% Input arguments : 
%
%   'eem' is a matrix of EEM data, with excitation wavelengths in the first
%   row, emission wavelengths in the first column, and corresponding
%   fluorescence intensity values for each ex/em wavelength combination.
%   The top left cell is ignored but copied to the output matrix in case it
%   contains information. Note that single emission and single excitation
%   scans are supported, but they must be in 'eem' format as described
%   (i.e. must include the excitation or emission wavelenth in the first
%   row or column, resp.)
%
%   'tol' is an n row by 2 colum matrix of correction tolerances in nm.
%   Columns specify tolerances below (left col) and above (right col) the
%   peak emission at each excitation wavelength.  Typical row assignements
%   are as follows :
%       row 1 : primary Rayleigh scatter
%       row 2 : primary Raman scatter
%       row 3 : secondary Rayleigh scatter
%       row 4 : secondary Raman scatter
%   Use tolerances >= 1000 to truncate corresponding portions of an EEM ad
%   set truncated values to 'baseopt' (note: truncation override
%   scising)
%   Typical example with truncation of emissions below excitation
%   wavelength:
%   [1000 12; ...
%    16   16; ...
%    18   18
%    18   18]
%
%   'coeff' is a matrix of polynomial curve fit coefficients which describe
%   the excitation wavelength-dependent emission wavelength for eah scatter
%   peak. Coefficients must be in rows and orderedd appropriately for the
%   'polyval' function (i.e. x^n, x^(n-1), ..., constant). Row assignements
%   must match the 'tol' matrix, but any power polynomial can be
%   speicified. Typical example for a second order polynomial matching
%   'tol' above:
%       [0        1.0000    0; ...
%        0.00058  0.8711    18.7770; ...
%        0        2.0000    0; ...
%        -0.0001  2.4085    -47.2964]
%   equivalent to: 0x^2 + 1x + 0, 0.0.00058x^2 + 0.08711x + 18.77, etc.
%
%   'baseopt' specifies optional baseline value to subsitute for values in
%   truncated scan regions. Default values is 0.
%
% Output parameters:
%
%   'eem_cor' is a matrix the dimensions as 'eem' containing corrected
%   fluorescence intensity values.
%   
%   'correct' is the correction matrix reflecting any changes following
%   validation
%   
%   'eem_filter' is the filter matrix used to correct 'eem' (i.e. a matrix
%   of ones and zeros the same dimensions as 'eem_cor', with zeros
%   representing values in scatter regions which were excised or truncated

%initialize outputs
eem_cor = [];
eem_filter = [];

if nargin >= 5 %check for minimum number of arguments
    if ~exist('baseopt') %assign default baseopt
        baseopt = 0;
    end
    
    if size(correct, 2) < 2 %single column format, replicate values
        correct = [correct correct];
    end
    
    cancel = 0;
    
    %validate input arguments
    %[ex, em, fl, fl_id] = unwrapeem(eem);
    fl = eem;
    em = em.';
    
    %ex = [ex(1) - (ex(2) - ex(1)), ex];
    %em = [em(1) - (em(2) - em(1)), em];
    
    if isempty(fl) %validate data matrix
        cancel = 1;
        errormsg = 'EEM data matrix is invalid';
    elseif size(parms, 1) ~= size(correct, 1) %validate polyfit parms
        cancel = 1;
        errormsg = 'Curve-fit coefficient matrix does not match correction matrix'
    end
    
    if cancel == 0 %proceed with analysis
        
        %form observed emission matrix
        em_obs = repmat(em, 1, length(ex));
        %assignin('base','em_obs',em_obs);
        %initialize filter matrices
        filt_excise = ones(size(em_obs));
        filt_trunc = zeros(size(em_obs));
        
        %loop through correction parameters for each peak
        for n = 1 : size(correct, 1)
            if correct(n, 1) > 0 | correct(n, 2) > 0 %test for no-filter condition
                peaks = polyval(parms(n, :), ex); %get array of scatter peak emissions
                
                pp=size(peaks)
                obb=size(em_obs)
                truc=size(repmat(peaks - correct(n, 1), length(em), 1))
                %form appropriate filter matrices for emissions below
                %scatter peak
                if correct(n, 1) < 1000 %excise
                    em_below = (em_obs - repmat(peaks - correct(n, 1), length(em), 1)) <= 0;
                else %truncate
                    %create logical index the same dimensions as em_obs
                    %with vals below the scatter peak lower limit = 1,
                    %update truncation filter
                    %[mobs]=size(em_obs)
                    %[m]=size(em)
                    
                    em_below = (em_obs - repmat(peaks, length(em), 1)) <= 0;
                    filt_trunc(em_below) = 1;
                end
                
                %form appropriate filter matrices for emissions above
                %scatter peak
                if correct(n, 2) < 1000 %excise
                    em_above = (em_obs - repmat(peaks + correct(n, 2), length(em), 1)) >= 0;
                else %truncate
                    %create logical index the same dimensions as em_obs
                    %with vals above the scatter peak lower limit = 1,
                    %update truncation fileter
                    em_above = (em_obs - repmat(peaks, length(em), 1)) >= 0;
                    filt_trunc(em_above) = 1;
                end
                
                %update excise filter matrix using combination of logical
                %indices (excise region = 0)
                filt_excise = filt_excise .* (em_below + em_above);
            end
        end
        
        %substitute 'baseopt' for truncated values
        fl(filt_trunc == 1) = baseopt;
        
        %update master filter to account for truncations
        filt_excise(filt_trunc == 1) = 1;
        
        %get index of excised values
        I_excise = find(filt_excise == 0);
        
        %excise values and grid data to reconstitute EEM if necessary
        if ~isempty(I_excise)
            %replace zeros with NaN
            fl(I_excise) = NaN;
            
            %form matched ex, em, fl vectors
            em_vec = reshape(em * ones(1, length(ex)), length(em) * length(ex), 1);
            ex_vec = reshape(ones(length(em), 1) * ex, length(em) * length(ex), 1);
            fl_vec = reshape(fl, length(em) * length(ex), 1);
            
            %get index of valid data points
            I_valid = find(~isnan(fl_vec));
            
            %interpolate using 2D or 3D algorithm as appropriate
            if size(fl, 2) > 1
                if size(fl, 1) > 1 %EEM matrix
                    %grid valid data points to form new EEM
                    fl = griddata(ex_vec(I_valid), em_vec(I_valid), fl_vec(I_valid), ex, em);
                else %excitation scan
                    %interpolate single ex scan
                    fl = interp1(ex_vec(I_valid), fl_vec(I_valid), ex, 'spline');
                end
            else %emission scan
                %interpolate single em scan
                fl = interp1(em_vec(I_valid), fl_vec(I_valid), em, 'spline');
            end
            
            %zero out any nulls or negative or near-zero values
            fl = nonneg(fl);
        end
        
        %assemble output EEM
        %eem_cor = wrapeem(ex, em, fl, fl_id);
        eem_cor = fl;
        
        %assemble total correction filter in EEM form
        filt_excise = filt_excise .* (filt_trunc < 1);
        eem_filter = wrapeem(ex, em, filt_excise);
        
    end
    
else
    errormsg = 'Too few arguments for function';
    cancel = 1;
end

if cancel == 1
    %clc
    disp(' '); disp(' ');
    disp(errormsg);
    disp(' '); disp(' ');
end

% function [ex, em, fl, fl_id, errormsg] = unwrapeem(eemdata)
% 
% % Syntax: [ex, em, fl, fl_id, errormsg] = unwrapeem(eemdata)
% %
% % Parses a consolidated matrix of eEM data and returns vectors of
% % excitation wavelengths ('ex') ans emission wavelengths ('em'), and a
% % matrix of fluroescence intensity values ('fl'). Any data in eemdata(1, 1)
% % is returned as 'fl_id'.
% 
% size(eemdata)
% 
% ex = [];
% em = [];
% fl = [];
% fl_id = [];
% errormsg = '';
% cancel = 0;
%                     
% if nargin > 0
%     %validate input
%     if size(eemdata, 1) <  2 | size(eemdata, 2) < 2
%         cancel = 1;
%         errormsg = 'Input matrix is invalid';
%         disp('A')
%     else
%         ex = eemdata(1, 2 : size(eemdata, 2));
%         em = eemdata(2 : size(eemdata, 1), 1);
%         fl = eemdata(2 : size(eemdata, 1), 2 : size(eemdata, 2));
%         fl_id = eemdata(1, 1);
% 
%         ex_inc = ex(2 : length(ex)) - ex(1 : length(ex) - 1);
%         I_ex = find(ex_inc <= 0);
% 
%         em_inc = em(2 : length(em)) - em(1 : length(em) - 1);
%         I_em = find(em_inc <= 0);
% 
%         if ~isempty(I_ex) | ~isempty(I_em) % check for nonmonotonic vectors
%             cancel = 1;
%             errormsg = 'Input matrix is invalid';
% 
%             ex = [];
%             em = [];
%             fl = [];
%             fl_id = [];
%         end
%     end
% else
%     cancel = 1;
%     errormsg = 'Too few arguments for function';
% end
% % 
% % 
function [eemdata, errormsg] = wrapeem(ex, em, fl, fl_id)
%syntax: [eemdata, errormsg] = wrapeem(ex, em, fl, fl_id)
%
% Concatenates vectors of excitation and emission wavelengths and
% corresponding fluorescence intensity matrix to produce a self-contains
% matrix of EEM data

eemdata = [];
errormsg = '';
cancel = 0;

if nargin >= 3
    ex_size = size(ex);
    em_size = size(em);
    fl_size = size(fl);
    
    %validate input
    if min(ex_size) > 1 | min(em_size) > 1 %confirm vector format
        cancel = 1;
        errormsg = '''ex'' and ''em'' must be vectors';
    elseif prod(fl_size) ~= (prod(ex_size) .* prod(em_size)) %test for mismatch
        cancel = 1;
        errormsg = 'Matrix dimensions do not match';
    else
        %test for nulls  in ex / em
        if length(find(isnan(ex))) > 0 | length(find(isnan(em))) > 0
            I_ex = find(~isnan(ex));
            ex = ex(I_ex);
            
            I_em = find(~isnan(em));
            em = ex(I_em);
            
            fl = fl(I_em, I_ex);
        end
        
        %check for non monotonic vectors
        ex_inc = ex(2 : length(ex)) - ex(1 : length(ex) - 1);
        I_ex = find(ex_inc <= 0);
        
        em_inc = em(2 : length(em)) - em(1 : length(em) - 1);
        I_em = find(em_inc <= 0);
        
        if ~isempty(I_ex) | ~isempty(I_em)
            cancel = 1;
            errormsg = '''ex'' and ''em'' must be monotonically increasing vectors';
        else %test for vector orientation
            if ex_size(1) > ex_size(2)
                ex = ex';
            end
            
            if em_size(2) > em_size(1)
                em = em';
            end
        end
    end

    if cancel == 0
        %use default ID label non specified
        if exist('fl_id') ~= 1
            fl_id = NaN;
        end
        
        %initialize zero matrix
        eemdata = zeros(length(em) + 1, length(ex) + 1);
        
        %assign values
        eemdata(1, 1) = fl_id;
        eemdata(1, 2 : size(eemdata, 2)) = ex;
        eemdata(2 : size(eemdata, 1), 1) = em;
        eemdata(2 : size(eemdata, 1), 2 : size(eemdata, 2)) = fl;
        
    end
else
    cancel = 1;
    errormsg = 'Too few arguments for function';
end

function m = nonneg(m)
% Converts matrix elements with values < 0 or NaN to 0

I = find(isnan(m));
m(I) = zeros(1, length(I)) + 1;
I = find(m <= eps);
m(I) = zeros(1, length(I)) + 1;
        