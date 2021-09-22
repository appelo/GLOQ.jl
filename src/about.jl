"""
    About GLOQ.jl
"""
function about()
    println("    ...`.`.`.`````.`````...`.`..``...`.`........`....................................................---
    `.````````.```.`````````.`..``.```.......`................`.......................................--
    .`````.`.``..```````````````..`............``............`.......................................---
    `````````````````````````.`.........----..........................................................--
    ````````````````````````...-/++oo+oyhhhhyyyssooosososssoso++osss+:................................--
    ``````````````````````...-/shddddddddddddddddddmmdmdddddddddddddddy+........`....``...............--
    ``````````````````````.:syhddmmdmmmmmmmmmmmmmmmmmmmmmmmmdmmmddmddmdds......`.```.``...............--
    ````````````````````.-ohhhddddmmdmmmmmmmmmmmmmmmmmmmmmmdddmddddddmmddo-.....``.`.```...............-
    ```````````````````./sshdddmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdho/:-....`````````...........-
    ``````````````````.+yyhddmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdys+:..``````.````..........
    `````````````````./hhhhddmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmdddho/-.``````.```.`.......
    `````````````````:hhhhhddmdmmdddddmmmmmddddddmmmmmmmmmmmmmdmdmdddddddmmmmddmmdds:..``````````.......
    `````````````````ohhhhhddddddhhhhhdhhhhhhhhddddddddddddddhhhhhhhhhhhhhddmmddmmmmdo-.``````````......
    ````````````````-yhhhhhhddddysssssyssysssyyyyyyyyyyyyyyyyyyyyyyyyyyyyyhddmmdmmmmmmdy-``````````.....
    ``````````````./yhhhhhhhddhysoo++ooo+++o++oooosoooooooossssssssssssyyyyhdmmmdmmmmmmmy.``````````.`..
    `````````````-sdddhhhhdddhyo+++++++++++/+++++++++++++++o+ooooooossssyyyyhdmmmmmmmmmmd/```````````...
    ````````````-yddddhhhddddyo+++++++++//////////////+++/+++++++ooooossssyyyhdmmdmmmmmmho````````````..
    ````````````ydddddhhhdddhs+++++++++////////////////++/+++++++++oooossssyyhhddddmmmmmdo`````````````.
    ```````````.hdddddhhhhhhyo+++++//+//////////////////////++++++++ooosssssyyhhddhdmmmmd/````````````..
    ```````````-ddddddhhyyyys++++++////////////////////+////++++++++ooosssssyyhhhdhdmmmmd/`````````````.
    ` `````````.yddddddhhyyso++++++/////////////////////////+//+++++ooosssssyyhhhhhhdmmmdo`````````````.
    ````````````-sddddddhyyso++++++///////////////////////////+/++++oooosssyyyyyyhdhdmmmd/``````````````
    ````` ```````/hhdhhhyysoo++++++////////////////////+//////+++++oooooossyyyyyyhhhddmdd:``````````````
    ```` ````````/hhhhhyysso+////+/+///////////////////++////++++ooooooosssyyyyyyyhhddddy.``````````````
      ` ```` ````-hhyhysso++/++//+++/+//////////////////++//+++++oosssssyhhhhhyyyyyyhhhh+.``````````````
        ` ````` `.shhyssoo+++/++++++ooooooo++o+++++++++++o+++ooossyyyhhhdddddhhyyyyyyhhhy/``````````````
      ` ` ` `` `./+syyosso+++//++++ossssyyyyyyyyysssysooosssyyyyyyyhhhhhhyhhhhhyhyyyyhyhy.``````````````
         `  ` ```/o++osssso++//++oooooossyyyyyyyyyyyyysoosyhhhhhhhhyyyyyyyyyhhhyyyhysyyh/```````````````
          ` `  ` `+o+ooooooo+++sooosssssssyyhhhhyyyyoysooosyhhhhhhhhdmmdhddhhhhyyyhyyyyy`````` ` ```````
               ` `:+++o+oo+++ssyosssyyhyyhddhyyyysossyo//++shhyyyyssyyyyyyyhhyyssyyyyys:`` `````` ``````
              `  ``/++++o++///+s+oossoooooooossssooshy+++++oydhyyysssosssyyyyssssshyyyo`` ` ``  ``` ````
                   :+ooo++++///s//++++++++oooosso++yys++/++oshhyysssssssssssosssssyyys/`   ``    `` ````
                  `-+oo++++++//+o///++++++++++++++oso+++/++osyhhysoooooooooossssssysss-`  `        `````
                   `/++++o+++///+o+++///////////+os+++++++++osyhhhsoooooooosssysssyys+`           ` ` ``
                   `.///oo+++++///ooo+///:////++sso+++++++++ooyhhyssssssssyysssssyyss+`      ``      ```
                     ://+oo++++//////+oooooo+++ooo++++++///+ooyyyyysooooooooossssyyss+        `     ` ``
                     .//+o++++++/////////////+ooo++++oo++++oosyhyyyysoooooooosssyyyss+`             ` ``
                     `/+oo+++++++//////////++o+///+ossssoossyhhdhhyyyssooooosssyyyyss/`              ` `
                      ./oo+++++++/+/////+++o++////++oooossyyhhhhhyyyssssoooosssyyyyss-`              ` `
                      `-+oo+++++++++///+++++//+++++++++ooosssssyyyyyyssssssssssyyyyso`               ```
                       ``/oo+o+++++++++++++++++++++++++ooooooossssyyyyyysooosssyyyysy:.`              ``
                    ``.-+sooooo++++++++///+++ooo+++++++++++ooosssyyhhyssoooossyyyyyhmmhs:.`          ` `
              ```.-/oyhdddhsooo+++++++///+++oooosssssssssssyyyyyhhhyysoo++oosyyyhhhmmmmmmhs:.``       `
       `  ``.:+syhdddddmmmmhsooo+++++++//++////+////++++++oooooossssoooo+oosyyyhhhhmNNNmmmmmds+-````  ``
    ```.:+syhddddddddmmmmmmmhssooo++++++/+///////++++ooooooosssssssoossooosyyyhhhhhhNNNNNNNmmmmdyo:.````
    +syddddddddddmmmmmmmmmmhysssoo+++++++++///////+++ooooosooooooooosssssssyhhhhhhhhdNNNNNNNNmmmmmmhs+-`
    dddddddddddmdmmmmmmmmmhssoosssoooo++++++////////+++++oooo+oooooossssyyyyhhdhhhhhhmNNNNNNNNNmmmmmdddy
    dddddmmmmmmmmmmmmmmmmyssooosssoooooo++++++///////+++++++ooooooosssyyyhhhhdddhhhddmNNNNNNNNNNNNmmmmmm
    mdmmmmmmmmmmmmmmmmmdsssoooosysssooooooo+++++/+////++++++ooosoosyyyyhhhhddhdhhhddddNNNNNNNNNNNNNNNNNN
    mdmmmmmmmmmmmmmmmmdsssooooosyyssssssssoooo++++++++oooooooossssyyhhhdhdddhddhhdddddNNNNNNNNNNNNNNNNNN
    mmmmmmmmmmmmmmmmmhsoooooooosyyosyyysssssssoooooooosssosossyyyyhhhdddddddhdhhddddddNNNNNNNNNNNNNNNNNN
    mmmmmmmmmmmmmmmdyooooooo++osyysosyyyyyyyyyyyyyssyyysssyyyhhhhhdddddddddhhhhhddddddmNNNMNMNNNMNNNNMMM
    mmmmmmmmmmmmmmdsoooooooo+ooyyyyoossyyyyyyyyyyyhyhhhhhhhhhhddhdhdhhhddddhhhhhdddddddNNMNNMMNNMNNMNMMM
    mmmmmmmmmmmmdhsooo+oooo+++osyyyoooosyyyyyyyyyyyyyhhhhhhhhhhhhdhhhhdddddhhhhdddddddhNNMMNMNNNMNNMMMMM
    mmmmmmmmmmmhsooooooooo++++osssss++oossyyyyyhyhhhhhhhhhhhhhhhhhhhhdddddhhhhhhddddhhhmNNMNMMMNMMNMMMMM
    mmmmmmmmmdyoooooooo+++++++ossysyo++ooosssyyyyhhhhhhhhhhhhhhhhhhhdddddhhhhhhhhhhhhyyhNNNMMMMMMMMNMMMM
    mmmmmmmdhsooooooo+oo++++++ossysss++++ooosssyyyyyyyyhhhhhhhhyyhhddddhhhhhhhhhhhyyyyyymNMNMNMNMNNMMMMM
    mmmmmdyssoooooo+++++++++++osssssyso++++ooooossssssyyyyyyyyyyyhdhhhhyyyhhyyyyyyyyyyyyhNMNMNMNMMNNMMMM
    mmmdysoooooooooooooo+++++oosssssyhso+++oooooossssssssssssyyyyhysyyysyyhyyyyyyyyyyssssmNNNMNNNMNMNMMM
    mmmdhhhyyssssoooooooooooosysossssyhso++oooooooooooooooosssyhyssssysssyyssssssysssssssyNNMMMNNNMNMMNM
    mmmmmmmmmmmmdddhhhhhhhddmmmhssssssyys+++++ooooooooo+ooossyhysooossooyyysssssssssssoooodNNNNNMNNMNNMM
    mmmmmmmmmmmmmmmmmmmmmmmmmmddossosssyyo++++ooo++++++++oosyysooooosooosyssosoooooooooooosNNNNNNNNNNMMM
    mmmmmmmmmmmmmmmmmmmmmmmmmmhyoosssssysso+++++++++////++oossooooosdhysssssooooooooooooooosNNMMNNNNNNMM
    mmmmmmmmmmmmmmmmmmmmmmmmmmhysssosssyssso++///////////+oossysyhdmmmNNmddhhyyyyhhhhhhhyyssdMNMNNNMNMMM")
    return nothing
end
