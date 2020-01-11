function [ctp, cycled_weights] = pcycler(Exp, Opt)
% [ctp] = pcycler(Exp, Opt)
% [ctp, cycled_weights] = pcycler(Exp, Opt)
%
% n.b. FIDs are untested as we don't use them in conventional pulse EPR.
%
% Calculates all possible coherence transfer pathways produced by a pulse
% sequence. This complete list is then pruned to remove all pathways which
% do not produce echoes (and optionally FIDs). These pathways and the 
% times at which they cross the desired pathway is calculated and
% outputted as ctp.
%
% If two outputs are requested, the weights of each echo after phase 
% cycling is also outputted (if no phase cycle is inputted all pulses 
% are assumed to have phase=0). The phase of the desired pathway to
% determine the phase of the detection. The real and imaginary components
% of all undesired pathways is calculated for each step in the phase cycle.
% This is outputted as cycled_weights, the rows of which correspond to ctp.
%
% Inputs
% Exp.nPulses       - Required, number of pulses to apply
% Exp.p_des         - Required, array length of Exp.nPulses + 1
%                     The desired coherence transfer pathway.
%                     n.b. first element must be rho0.
% Exp.d_names       - Required, symbollic array length Exp.nPulses - 1
%                     The names of the time delays.
% Exp.var           - Required, symbol. Name of incremented time delay.
%                     Used to calculated relative times of echo crossings.
% Exp.p_cycle       - Optional, n_step x Exp.nPulses array specifying
%                     phases of pulses applied in RADIANS !!
% Exp.max_order     - Optional, maximum coherence order to examine
%                     default Exp.max_order = 1.
% Exp.p0            - Optional, coherence order at start of sequence
%                     default Exp.p0 = 0.
% Opt.FID           - Optional, if Opt.FID = 0 (default) FIDs will not be
%                     considered as valid pathways.
% Opt.onlycrossing  - Optional, if Opt.noncrossing = 0 (default) echos
%                     which do not cross the desired pathway are  output
%                     regardless. If set to 1, only echos which cross are
%                     outputted.
% Opt.filter        - Optional, array length Exp.nPulses
%                     Used for coherence filtering. To allow any pathway
%                     set to 1i*ones(1, Exp.nPulses). Each
%                     element corresponds to a delay. If any element of
%                     Opt.filter is real the filter is activated, and 
%                     coherence orders not equal to those specified are
%                     ignored.
%                     e.g. To remove all coherences during the third delay
%                          of a 5 pulse sequence (as in RIDME) set:
%                              Opt.filter = [1i, 1i, 0, 1i, 1i]
%                     Default is 1i*ones(size(Exp.nPulses)).
%                     Note that the requirement of coherence order -1 for
%                     detection is effectively equivalent to a filter
%                     where the last element is -1.
%
% Outputs
% ctp               - cell array, nEchos x nPulses + 2
%                     Coherence transfer pathways (starting at p0) which
%                     lead to echos. The final column is a character array
%                     with the time at which the echo crosses the desired
%                     coherence transfer pathway. If the echo does not
%                     cross then it is denoted n[i] where i is the number
%                     of the echo. See Opt.onlycrossing.
% cycled_weights    - complex array, length nEchos
%                     The real and imaginary parts of each coherence
%                     transer pathway after phase cycling.

% Set defaults and ensure input correct
% nPulses
if ~isfield(Exp, 'nPulses')
    error("Field Exp.nPulses is not optional")
end
if ~isscalar(Exp.nPulses) || floor(Exp.nPulses) ~= Exp.nPulses
    error("Field Exp.nPulses must be a scalar integer")
end

% Maximum coherence order
if ~isfield(Exp, 'max_order')
    Exp.max_order = 1;
else
    if ~isscalar(Exp.max_order) || floor(Exp.max_order) ~= Exp.max_order
        error("Field Exp.max_order must be a scalar integer")
    end
    % Ensure positive
    if Exp.max_order < 0
        Exp.max_order = -Exp.max_order;
    end
end

% Initial coherence order
if ~isfield(Exp, 'p0')
    Exp.p0 = 0;
else
    if ~isscalar(Exp.p0) || floor(Exp.p0) ~= Exp.p0
        error("Field Exp.p0 must be a real scalar")
    end
    % Ensure initial coherence order is possible
    if abs(Exp.p0) > Exp.max_order
        error("Initial coherence order Exp.p0 greater than maximum " + ...
               "coherence order specified by Exp.max_order.")
    end
end

% Desired coherence pathway
if ~isfield(Exp, 'p_des')
    error("Field Exp.p_des is not optional")
end
if floor(Exp.p_des) ~= Exp.p_des
    error("Exp.p_des must be an integer array")
end
if ~isvector(Exp.p_des)
    error("Exp.p_des must be a vector array")
end
if Exp.p_des(1) ~= Exp.p0
    error("Desired coherence transfer pathway must start with Exp.p0")
end
if Exp.p_des(end) ~= -1
    error("Desired coherence transfer pathway must end with -1")
end
% Check that all specified orders are possible
if sum(abs(Exp.p_des) > Exp.max_order) ~= 0
    error("No element of Exp.p_des can be larger than Exp.max_order")
end
% Make sure it's the right length and is row vector
if size(Exp.p_des, 1) ~= 1
    Exp.p_des = Exp.p_des.';
end
if size(Exp.p_des, 2) ~= Exp.nPulses + 1
    error("Exp.p_des must be Exp.nPulses + 1 in length")
end

% Names of delays
if ~isfield(Exp, 'd_names')
    error("Field Exp.d_names is not optional")
end
if ~isa(Exp.d_names, 'sym') % CHECK THIS WORKS
    error("Exp.d_names must be an array of symbols")
end
if ~isvector(Exp.d_names)
    error("Exp.d_names must be a vector array")
end
% Make sure it's the right length and is row vector
if size(Exp.d_names, 1) ~= 1
    Exp.d_names = Exp.d_names.';
end
if size(Exp.d_names, 2) ~= Exp.nPulses - 1
    error("Exp.d_names must be Exp.nPulses - 1 in length")
end

% Phase cycle - only checks size
if ~isfield(Exp, 'p_cycle')
    Exp.p_cycle = zeros(1, Exp.nPulses);
end
if size(Exp.p_cycle, 2) ~= Exp.nPulses
    error("Each row of Exp.p_cycle must be Exp.nPulses in length")
end

% If Opt not inputted set all to defaults
if nargin == 1
    Opt.FID = 0;
    Opt.onlycrossing = 0;
    Opt.filter = 1i*ones(1, Exp.nPulses);
else
    % Include FIDs ?
    if ~isfield(Opt, 'FID')
    	Opt.FID = 0;
    end
    if Opt.FID ~= 0 && Opt.FID ~= 1
        error("Opt.FID must be 0 or 1.")
    end
    % Ignore echos which do not cross?
    if ~isfield(Opt, 'onlycrossing')
    	Opt.onlycrossing = 0;
    end
    if Opt.onlycrossing ~= 0 && Opt.onlycrossing ~= 1
        error("Opt.onlycrossing must be 0 or 1.")
    end
    % Filtering
    if ~isfield(Opt, 'filter')
    	Opt.filter = 1i*ones(1, Exp.nPulses);
    end
    if size(Opt.filter, 1) ~= 1
        Opt.filter = Opt.filter.';
    end
    if size(Opt.filter, 2) ~= Exp.nPulses
        error("Opt.filter must be a vector length Exp.nPulses")
    end
end

% Start doing maths
% Generate all valid coherence transfer pathways
p_t = Exp.p0;
for p = 1:Exp.nPulses
    % Apply pulse
    p_t = apply_pulse(p_t, Exp.max_order);
end

% Remove all which do not end in -1 to avoid duplication
cc = p_t(:, end) == -1;
p_t = p_t(cc, :);

% Interpulse delays
delays = sym('t', [1, Exp.nPulses]);    
% Blank coherence times
t_echo = sym(-1i.*ones(length(p_t), 1));
t_FID = t_echo;
cc = zeros(length(p_t), 1);
dd = cc;

% Find all echos and FIDs
for p = 1:length(t_echo)
    
    % Echos spend equal time in the +1 and -1 coherence orders
    % Number of periods of evolution with p = +1
    n_pos = sum(p_t(p, :) == +1);
    % Number of periods of evolution with p = -1
    n_neg = sum(p_t(p, :) == -1);
    
    % If an echo is possible
    if (n_pos >= 1) && (n_neg >= 1)
        % Remove placeholder
        t_echo(p) = 0;
        cc(p) = 1;
        % Loop over pulses
        for i = 1:Exp.nPulses
            % Create echo time
            if p_t(p, i+1) == 1
                t_echo(p) = t_echo(p) + delays(i);
            elseif  p_t(p, i+1) == -1
                t_echo(p) = t_echo(p) - delays(i);
            end
        end
    end
    
    % FIDs occur when magnetisation changes from 0 to -1.
    % (FIDs with coherence order +1 are ignored).
    % The FID component of pathways which can give echos is discounted.
    if Opt.FID == 1 && cc(p) == 0
        n_FID = 0;
        for q = 1:Exp.nPulses
            if p_t(p, q) == 0 && p_t(p, q + 1)
                n_FID = n_FID + 1;
            end
        end
        % If only one FID (more than one is considered in echos section)
        if n_FID == 1
            t_FID(p) = sum(delays(p_t(p, 2:end) ~= -1));
            dd(p) = 1;
        end
    end
end

% Remove pathways which do not lead to echoes (keep FIDs if required)
if Opt.FID == 0
    t_coh = t_echo(cc == 1);
    p_coh = p_t(cc == 1, :);
else
    t_coh = [t_echo(cc == 1); t_FID(dd == 1)];
    p_coh = [p_t(cc == 1, :); p_t(dd == 1, :)];
end

% Apply coherence filter if applied
for d = 1:Exp.nPulses
   if isreal(Opt.filter(d))
       cc = p_coh(:, d+1) == Opt.filter(d);
       t_coh = t_coh(cc);
       p_coh = p_coh(cc, :);
   end
end

% Find desired pathway
coh_ind = 0;
for e = 1:length(t_coh)
    if p_coh(e, :) == Exp.p_des
        coh_ind = e;
        break
    end
end
if coh_ind == 0
    error("Desired coherence pathway could not be found in list of" + ...
          "possible pathways, please check p_des.")
end

% Observed echo time = time up to last pulse + acquisition delay
% FIDs do not need this as do not depend on delay time
sym temp;
for e = 1:length(t_coh)
    temp = solve(t_coh(e) == 0, delays(end));
    if size(temp) ~= 0
        t_coh(e) = sum(delays(1:end-1)) + temp;
    end
end

% Convert symbols to conventional notation
for d = 1:length(delays) - 1
    t_coh = subs(t_coh, delays(d), Exp.d_names(d));
end

% Find location of detection
coh_loc = t_coh(coh_ind, :);

% Find crossing time by solving for time variable
t_coh = t_coh - coh_loc;
sym temp;
t_cross = sym('n', size(t_coh));
for e = 1:length(t_coh)
    temp = solve(t_coh(e) == 0, Exp.var);
    if size(temp) ~= 0
        t_cross(e) = temp;
    end
end

% Readable output
ctp = num2cell(p_coh);
for i = 1:size(p_coh, 1)
    ctp(i, size(p_coh, 2) + 1) = {char(t_cross(i))};
end

% If requested, apply phase cycle
if nargout == 2
    signal_phase = zeros(size(p_coh, 1), size(Exp.p_cycle, 1));
    dp = zeros(1, Exp.nPulses);

    % Loop over coherence transfer pathways
    for e = 1:size(p_coh, 1)

        % Get delta p
        for p = 1:Exp.nPulses
            dp(p) = p_coh(e, p + 1) - p_coh(e, p);
        end

        % Loop through phase cycles
        for c = 1:size(Exp.p_cycle, 1)
            % Start at -y
            signal_phase(e, c) = -pi;
            % Loop through steps
            for p = 1:Exp.nPulses
                % Calculate phase shift
                signal_phase(e, c) = signal_phase(e, c) - dp(p)*Exp.p_cycle(c, p);
            end
        end
    end


    % Cancel unwanted echos
    cycled_weights = zeros(1, size(p_coh, 1));
    for e = 1:size(p_coh, 1)
        for c = 1:size(Exp.p_cycle, 1)
            phase = exp(1i*(signal_phase(e, c) - signal_phase(coh_ind, c)));
            cycled_weights(e) = cycled_weights(e) + phase;
        end
    end
    
    % Readable output (chop out small values and normalise)
    cycled_weights = cycled_weights.';
    for e = 1:size(cycled_weights)
        re = real(cycled_weights(e));
        im = imag(cycled_weights(e));
        if abs(re) <= 1e-15; re = 0; end
        if abs(im) <= 1e-15; im = 0; end
        cycled_weights(e) = re + 1i*im;
    end
end

end


function order_out = apply_pulse(order_in, max_order)
% Applies a pulse to an array of coherence orders. The final column of
% order_in is taken to be the coherence orders before the pulse. order_out
% contains all possible coherence order pathways after the pulse.

% Number of orders to loop over
n_in = size(order_in, 1);
% Number of coherence orders up to and including current
n_p = size(order_in, 2) + 1;

% Output array with imaginary in first column
% Anything left imaginary after loop is removed
order_out = zeros(5*n_in, n_p);
order_out(:, 1) = -1i;

% Set ticker
tic = 0;

% Loop over inputted coherences
for c = 1:n_in
    
    % Loop over possible dp values
    for dp = -2*max_order:1:2*max_order
        
        % New coherence order
        temp = order_in(c, end) + dp;
        
        % Move index
        tic = tic + 1;
        
        % If coherence order feasible, store initial and final values
        if abs(temp) <= max_order
            order_out(tic, 1:n_p - 1) = order_in(c, :);
            order_out(tic, end) = temp;
        end
        
    end
    
end

% Cut out unfeasible orders
keep = order_out(:, 1) ~= -1i;
order_out = order_out(keep, :);

end