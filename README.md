# PCycler
PCycler is a MATLAB function which aids the creation of efficient phase cycles in pulse magnetic resonance experiments. For an arbitrary pulse sequence, PCycler can apply a virtual phase cycle and outputs which echos and FIDs cross the desired resonance feature and at what times. This is useful for pulse Electron Paramagnetic Resonance (EPR) and Nuclear Magnetic Resonance (NMR).

## Setup
PCycler requires the [Matlab Symbolic Toolbox](https://uk.mathworks.com/products/symbolic.html).

To install PCycler, you need to add the file pcycler.m to your Matlab path, see [here](https://uk.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) for the official Matlab documentation on how to do this.

## Usage
A full list of inputs and outputs is detailed in the header of [pcycler.m](pcycler.m). An example input script is included ([DEER4P_2Step.m](DEER4P_2Step.m)).

## License
This project is licenced under the GNU General Public License, see [LICENSE](LICENSE) for more details.

## Authors
Edmund Little

## Acknowledgements
This function is a practical implementation of:

Phase Cycling in Electron Spin Echo Envelope Modulation, _S. Stoll and B. Kasumaj_, Appl. Magn. Reson. (__2008__) 35, 15-32

Copyright (c) 2019: Edmund Little
