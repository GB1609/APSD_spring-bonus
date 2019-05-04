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
    to_write="";
  }

  void update_string(string to_append)
  {
    to_write+=to_append;
  }

  void update_string(int nt,double ser, double par,double su,string type)
  {
    string to_append = to_string(nt) + ";" + to_string(ser) + ";" + to_string(par) + ";" + to_string(su) + ";" + type+"\n";
    to_write=to_write+to_append;
  }

  void write()
  {
    ofstream file_ex;
    file_ex.open ("sheets/"+name_file+".txt",ios::out |ios::app);
    file_ex << to_write;
    file_ex.close();
  }
  void write(string s)
  {
    ofstream file_ex;
    file_ex.open ("sheets/"+name_file+".txt",ios::out |ios::app);
    file_ex << s<<"\n";
    file_ex.close();
  }
};

#endif /* FILEWRITER_H_ */
