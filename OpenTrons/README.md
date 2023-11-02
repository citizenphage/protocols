Opentrons protocols used by the CPL for automating analyses.

Pipeline overview:

<img width="805" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/dec80adc-7fe6-4ba2-8fb6-515e7cad8854">

Our pipeline offers both a semi-solid and a liquid-based approach to phage isolation:
<img width="790" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/aa0b176f-016a-4beb-a957-e1dd94107d12">
<img width="795" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/7b5ea0d4-2298-4456-b375-898331ce45f3">

Each protocol folder contains the protocol template file, accompanied by by a usage guide.

These protocols contain values for the user to change, found near the top of the python script under the heading "USER DEFINED VARIABLES". These values are required for customisation of the protocol run. They include details such as:
    -Number of wells containing sample
    -Starting volume of media 
    -Optional parameters

Once these variables are defined, we suggest saving the protocol with a name that allows run tracibility. You may wish to include the date or daycode, the protocol ID, and the run number (protocol can be used more than once in a day). Example file name: "230821_CPL1_RUN1.py".


Robot usage:
The Opentrons OT-2 liquid handling robot contains 11 deck slots and a tip bin.
<img width="424" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/3fc713e7-10bf-40cd-b9a3-a582a0227860">

It is important to regularly calibrate the robot for it to move and interact with labware accurately.
Keep the robot door closed whenever possible. All internal surfaces should be disinfected after use.
Tips are used from top to bottom, left to right.

Also provided within this repository is a 3D-printable model of a 50ml tube rack. This is OT-2 compatible and the only tube rack required for this pipeline. When using, ensure the holes for the tubes are on the right hand side of the labware.
<img width="640" alt="image" src="https://github.com/citizenphage/protocols/assets/101196413/049297bd-1a95-4302-8805-0a0d472f025a">

