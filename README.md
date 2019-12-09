# mininap-N-gap-filling-tools-
N-gap filling tools using minimap


//  Created by Walfred MA in 2018, wangfei.ma@ucsf.edu.
//  Copyright © 2018 UCSF-Kwoklab. All rights reserved.


---------description---------

N-gap filling tools using blastn

1.it detect Ngaps or low quality regions in query sequences.
2.it uses 500 sizes as anchor to locate Ngaps on reference assembly.
3.it replace Ngaps or low quality regions with it anchors with sequences from reference assembly.

---------usage---------

Python2.7 -i querypath -r refpath -o outputpath
