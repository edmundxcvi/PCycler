# PCycler
Calculate residual coherences in magnetic resonance experiments

[ctp] = pcycler(Exp, Opt)
[ctp, cycled_weights] = pcycler(Exp, Opt)

n.b. FIDs are untested as we don't use them in conventional pulse EPR (for which this script was written).

Calculates all possible coherence transfer pathways (ctp) produced by a pulse sequence. This complete list is then pruned to remove all pathways which do not produce echoes (and optionally FIDs). These pathways and the times at which they cross the desired pathway is calculated and outputted as ctp.

If two outputs are requested, the weights of each echo after phase 
cycling is also outputted (if no phase cycle is inputted all pulses 
are assumed to have phase=0). The phase of the desired pathway to
determine the phase of the detection. The real and imaginary components
of all undesired pathways is calculated for each step in the phase cycle.
This is outputted as cycled_weights, the rows of which correspond to ctp.

Inputs
Exp.nPulses       - Required, number of pulses to apply
Exp.p_des         - Required, array length of Exp.nPulses + 1
                    The desired coherence transfer pathway.
                    n.b. first element must be rho0.
Exp.d_names       - Required, symbollic array length Exp.nPulses - 1
                    The names of the time delays.
Exp.var           - Required, symbol. Name of incremented time delay.
                    Used to calculated relative times of echo crossings.
Exp.p_cycle       - Optional, n_step x Exp.nPulses array specifying
                    phases of pulses applied in RADIANS !!
Exp.max_order     - Optional, maximum coherence order to examine
                    default Exp.max_order = 1.
Exp.p0            - Optional, coherence order at start of sequence
                    default Exp.p0 = 0.
Opt.FID           - Optional, if Opt.FID = 0 (default) FIDs will not be
                    considered as valid pathways.
Opt.onlycrossing  - Optional, if Opt.noncrossing = 0 (default) echos
                    which do not cross the desired pathway are  output
                    regardless. If set to 1, only echos which cross are
                    outputted.
Opt.filter        - Optional, array length Exp.nPulses
                    Used for coherence filtering. To allow any pathway
                    set to 1i*ones(1, Exp.nPulses). Each
                    element corresponds to a delay. If any element of
                    Opt.filter is real the filter is activated, and 
                    coherence orders not equal to those specified are
                    ignored.
                    e.g. To remove all coherences during the third delay
                         of a 5 pulse sequence (as in RIDME) set:
                             Opt.filter = [1i, 1i, 0, 1i, 1i]
                    Default is 1i*ones(size(Exp.nPulses)).
                    Note that the requirement of coherence order -1 for
                    detection is effectively equivalent to a filter
                    where the last element is -1.

Outputs
ctp               - cell array, nEchos x nPulses + 2
                    Coherence transfer pathways (starting at p0) which
                    lead to echos. The final column is a character array
                    with the time at which the echo crosses the desired
                    coherence transfer pathway. If the echo does not
                    cross then it is denoted n[i] where i is the number
                    of the echo. See Opt.onlycrossing.
cycled_weights    - complex array, length nEchos
                    The real and imaginary parts of each coherence
                    transer pathway after phase cycling.
