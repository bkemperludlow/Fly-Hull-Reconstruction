Hey Sam,

Here is the cine compression package (sparser).

I did not include the Phantom SDK. You need to load the Phantom DLLs first. 

We put them in "cine2sparse\PhantomSDK\SDK13.4.787.0"
and that's the version we use.

then to load Phantom Libraries from that location run "runMeFirst.m"
or use your own loader.

To compress cine files use Cine2Sparse.m (it can do batch conversion by selecting multiple files).

To visualize sparsed movies use SparseMovieLoaderV2.m
It can visualized three files at a time, just select multiple file after clicking the "Load" button on the bottom left. It has other buttons to play with as well (note the little arrow buttons on the two sides of the slider).

This viewer requires the toolbox "GUI Layout Toolbox", which can be loaded from the Matlab Addons menu (get add-ons).

reading an image from the sparse structure is done using ImfromSp.m

Please let me know if you have any questions.

Good luck,
Tsevi


