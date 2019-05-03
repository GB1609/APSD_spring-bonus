/*
 * FileWriter.h
 *
 *  Created on: 3 mag 2019
 *      Author: gb1609
 */

#ifndef FILEWRITER_H_
#define FILEWRITER_H_

class FileWriter
{
  private:
  int excercise;
  string name_file;
  string to_write;
  public:
  FileWriter(int ex):excercise(ex)
  {
    name_file="excercise"+to_string(excercise);
    to_write=name_file+"\n";
  }

  void update_string(string to_append)
  {
    to_write+=to_append;
  }

  void write()
  {
    ofstream myfile;
    myfile.open (name_file+".txt",ios::out |ios::app);
    myfile << to_write;
    myfile.close();
  }
};

#endif /* FILEWRITER_H_ */
