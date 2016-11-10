//
// Created by syd on 16-10-23.
//

#ifndef LIBMBR_MBR_H
#define LIBMBR_MBR_H

int mbr_encode(char *data);
int mbr_decode(int len,char **data);
int mbr_regenerate(char *data);



#endif //LIBMBR_MBR_H
