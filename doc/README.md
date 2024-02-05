# Creating Documentation Using ROFF Language

To generate documentation using the ROFF language, 
it is necessary to create a dedicated file for each command. 
For instance, create.1. After writing the file, 
compile it using the `groff` software available in the Linux repositories. 
Execute the following command in a Linux 
console whenever the documentation for a particular command is updated:

```bash
groff -man -Tascii create.1 > create.man
