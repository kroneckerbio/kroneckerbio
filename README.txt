KroneckerBio is a modelling toolbox for systems biology written in Matlab. It is maintained by David R Hagen in the Bruce Tidor lab at MIT.

Quick start:
Run the InitKronecker.m script. This will add appropriate folders to your Matlab path. The easiest way to learn how to use the toolbox is to go to the Testing folder and try out the scripts. The functions are in the Source folder. Functions NamedLikeThis are top level functions (look here for useful stuff), those namedLikeThis are helper functions (mostly used by the top level functions), and those namedlikethis are mathematical function to supplement basic Matlab functionality (useful even outside Kroncker).

Help:
The help information on the functions is mostly up-to-date. Contact David Hagen (david@drhagen.com) for help, to report bugs, or to request features.

History:
KroneckerBio began as a test to see if mass action modelling could be done in a way to substantially improve performance and ease-of-use of systems biology tools. Performance is improved because purely mass action models can be expressed mathematically as simple matrix multiplications, which are very fast in Matlab. Ease-of-use is improved because mass action reactions are easy to write (A + B -> C) while still being able to express the mechanistic behavior of biological systems. Prof. Bruce Tidor and Prof. Jacob White were the PIs behind the project with Joshua Apgar and Jared Toettcher as first developers. David Hagen took over as developer and totally rewrote it.
