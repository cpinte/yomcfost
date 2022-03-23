/* Author : Guieu Sylvain.
 * $Log: openASCII_v2.i, v 2 $
 *
 * 2006/06/09 14:35  guieu
 * openASCII_v2.i added to contrib
 * 2006/09/06  guieu
 * Add the : include "string.i"
 * And change strjoin to openascii_strjoin
 *  
 * -----  v2.1 ------
 * Add possibility to read several table in the same file.
 * Add possibility to separete columns with a char  (sep keyword)
 *  
 */

#include "string.i"

if(_PRINT) write,"include \"openASCII_v2.i\"";


/*  ************* Extern Variable ******************** */
_default_ASCII_format_file = "~/yorick/open_ASCII_default_format.dat" ;

openascii_comment_keys = ["comment"];
openascii_comment_string ="#";
openascii_format_def_string = "!";
openascii_sep_char = ""; /* default separation char to read table,
                            "" is whitespace */


_max_alt_key_number = 100;
_max_key_in_header  = 5000;
_PROMPT =1;

/* **************** structure definition ************ */
struct KEYS {
  string name;
  string  type;
  int num;
  int pos;
}
struct ALT_NAME {
  string name; 
  int num;
  int pos;
}
struct OPERATION {
  string name;
  int    col1;
  int    col2;
  string op;
  int pos;
}

  func OpenASCII2 (file, &f, N=, prompt=,
                 file_format=, structure=,
                 add_kw=, debug=, code_file=, no_op=,
                  onlystruct=, noheader=, struct_name=, valtype=, tab=, sep=)
{
  /* DOCUMENT
   struct_out = OpenASCII (file)
                OpenASCII (file, f, tab=num)
   OpenASCII (version 2) 
   
   OpenASCII read the ASCII file FILE and put results in a structure.

   If there is several table in the same file (separated by at least one void
   line) You can access to the table NUM with the tab keyword.
   f is a optional output text stream ready to read a next table if exist.
   You can do :
      f = open("my_file");
      s1 = ReadASCII(f);  // read the first table
      s2 = ReadASCII(f);  // read the second table ...
      s4 = ReadASCII(f,tab=2); // skip the third table and read the forth
      
      
   The structure extentions are defined by the header of the file and an
   other special file wich possess definition of structure you want to built.
   specify the default path of this file with the   
   _default_ASCII_format_file    extern variable

   keyword N specify the number of data lines to read, if N exceed the table size,
   the structure returned have a dimension of N but element in sus are not
   filled and a warning message is printed.
   
   
   The file format is like this :
   
   extension_1 format_1 key_1_1 
   extension_2 format_2 key_2_1 key_2_2 key_2_3
   extension_3 format_3 key_3_1
      .          .           .           .          .
      .          .           .           .          .
      .          .           .           .          .
   extension_m format_m key_m_1 key_m_2 ... ... ... key_m_n

   Where :
      extension_i : the structure extension of the element i
                    must be a valid yorick varaible ( A-Z , a-z,  _ )
      format_i    : the format of the element i must be : 
                    string, double, float, int, long,
                    char or special (see bellow)
      key_i_n     : Key-words of the element. If there is one this key-words
                    in the file header the structure element is create and
                    filed with the corresponding column. key-words do not
                    contained any spaces. you can specify a number of
                    key-word you want. 

   A line like : " extention type"
                                               is equivalent to
                  "extention type extention"
                  key-words and extenction are the same string.
                  
   You can add lines with the ADD_KW keyword (vector of string of
   format definition)
   Or you can add line in the header of the file, this line start with the
   "!" string (default) you can change this string
   with the openascii_format_def_string extern variable.
   
   A line added like this have priority on existing extenction definition
   on the foramt_file. Une can also use this if you change temporaly one of
   varaible types.

   you can set  _default_ASCII_format_file=""  or
   OpenASCII("my_file", file_format="")
   If you want that OpenASCII do not use a file format but only definitions
   in the file (started with "!")
   
   If 2 different lines have the same key-word, the 2 extentions will be
   created in the structure, but the type (string, float, ... etc)
   used to read the file is the most on the left in the line definition.
   
   Exemple :
   _______________________________
   |name string name NAME id      |
   |num     int   id   ID  numero | 
   |______________________________|
   Here, for read 'id' in the header of file the second line win and id
   column (in file) is read  with  the 'int' type. 'name' extenction
   is created on the structure but will not be filled if is type is
   incompatible (the case here)
   But :
    _______________________________
   |name string   id    NAME name  |
   |num     int   id    ID   numero|  
   |_______________________________|
   Here 'id' have the same positions the first line win because it declared
   before. 

   IF THE (WINNER) TYPE DO NOT MATCH THE REAL COLUMN TYPE
   (EX DOUBLE INSTEAD OF STRING) THE READ FUNCTION FAILED.
   A WARNING message is printed. 

   
   Example of format file 
   ____________________________________________________
  | name string	 name	NAME	Noms                   |
  | Ra	 string	 ra	RA	Ra                     |
  | Dec	 string	 dec	DEC	Dec                    |
  | T	 double	 t	T	Teff	Temp           |
  | L	 double	 l	L	Lbol	lbol	L/Ls   |
  | v	 double	 v	V	mV	Mv             |
  | r	 double	 r	R	mR	Mr             |
  | i	 double	 i	I	mI	Mi             |
  | z	 double	 z	Z	mZ	Mz             |
  | j	 double	 j	J	mJ	Mj             |
  | h	 double	 h	H	mH	Mh             |
  | k	 double	 k	K	mK	Mk             |
  | rad  double  radus  rayon   R/Rsun  R              |
  | id   int     id     ID                             |
  |____________________________________________________|                

  The tabe to read is like this :
   _____________________________________________________________________
  |# some comment                                                       |
  |! pata float Pata                                                    |
  |# other comments                                                     |
  |#id  name             RA              tsp     T     Pata  comment    |
  |364  KPNO-Tau-1	04:15:14.714 	9.00	2571  1.2   brown dwarf |
  |370  KPNO-Tau-6	04:30:07.244 	9.00	2571  2.3               |
  |396  KPNO-Tau-4	04:27:27.997 	9.50	2500  4.6   very young  |
  |390  KPNO-Tau-10	04:17:49.554 	5.50	3065  0.3               |
  |_____________________________________________________________________|
  
  All the comment lines must be at the start of the file and must be started
  with the "#" string (default). You can change it
  with the openascii_comment_string extern variable.
  If there is a void line OpenASCII assume that there is an other table.

  THE LAST COMMENT LINE IS THE HEADER
  
  In this exemple the line "! pata float Pata" add the definition of
  "Pata". Because you do not use often the "Pata" key-word and you just
  want to define it for this file. Or,  "Pata" is defined in your
  format file like a string and you want defined it for this file like a float.
  
  If the LAST key-word in the header is "comment"  or any of
  string defined in the openascii_comment_keys extern N dimension variable,
  OpenASCII use a specific format to read this column. This column can
  contain void string or string with any character (" " also).
  The end of lines is read entirely. 
  
  You can add the definition of this kind of variable with the
  "special" format  example :
  -----------------------------------------------
  com special comment COMMENT COM END
  ------------------------------------------------
  If one of comment, COMMENT, COM or END is at the and of Header the column
  will be read specialy (see above), else this is a simple string variable. 

  Also you can separate EACH columns whith an other thing than whitespace
  using the SEP keyword.
  ------------------------------------
  #id	| name	      |	 tsp  | T     
  364	| KPNO Tau 1  |	 9.00 | 2571  
  370	| KPNO Tau 6  |	 9.00 | 2571  
  396	| KPNO Tau 4  |	 9.50 | 2500  
  390	| KPNO Tau 10 |	 5.50 | 3065  
  ____________________________________
  Can be read with   s = OpenASCII(my_file, sep="|" )
  
  
  OpenASCII look if some variable can be computed with other ones.
  example : if there is in the header : I  and I-Z, OpenASCII
  compute  Z = I - (I-Z)  and consider Z like a new column.  
  You can switch of this with the NO_OP keyword. 
  
  User can force the used structure with the STRUCTURE keyword.  only
  structure extension found by OpenASCII will be filed.  If you give a
  structure definition and ONLYSTRUCT keyword is non nul, OpenASCII use
  only the definition of structure to read the file, no definition
  will be read in your format file, the keys-word in header must
  be the same string of structure extension.
  If NOHEADER keyword is non nul or if there is no header in the file
  and if you have specified a structure definition, openASCII will try
  to read the file with the structure definition, order of columns will
  correspond to order of elements in structure definition. 

  
  
  

  OpenASCII create a temporally file with the structure definition
  and the function definition readable by human. By default OpenASCII
  include this file in yorick and erase it, you can specify the code
  file with the CODE_FILE keyword and change it by hand if there is some
  problems or if you want to create a template function for read several file
  with same cstructure. 
  The function created is named tmp_func_openASCII and is erased and replaced
  when OpenASCII is called again.
  If you do :
  > test = OpenASCII("my_file.dat");
  The tmp_func_openASCII is created and openASCII return the result of
  tmp_func_openASCII();
  So you can redo :
  
  KEYWORD :
  file_format , prompt ,  structure, add_kw , code_file
   SEE ALSO:
 */
  local f;
  f = open(file);
  return ReadASCII2 (f, N=N, prompt=prompt,
              file_format=file_format, structure=structure,
              add_kw=add_kw, debug=debug, code_file=code_file,
              no_op=no_op,
              onlystruct=onlystruct,
              noheader=noheader, struct_name=struct_name,
              valtype=valtype, tab=tab, sep=sep);

}
  
func ReadASCII2 (&f, N=, prompt=,
                 file_format=, structure=,
                 add_kw=, debug=, code_file=, no_op=,
                 onlystruct=, noheader=, struct_name=, valtype=, tab=, sep=)

{
  extern tmp_func_openASCII;
  /* erase the previoue tmp_func */
  tmp_func_openASCII = [];
  if (is_void(func_name)) func_name = "tmp_func_openASCII";
  
  if (is_void(prompt)) prompt = _PROMPT;
  if (onlystruct && is_void(structure))
    error, "You must specify a structure definition with onlystruct option";
  //code_file = file+".i";
  if (is_void(sep)) sep=openascii_sep_char;
  sep=strtrim(sep);
  sep = strpart(sep,1:1);
  
  if (!is_void(code_file)) debug=1;
  
  comment_string = (openascii_comment_string ? openascii_comment_string : "#");
  format_def_string=(openascii_format_def_string ? openascii_format_def_string : "!");
  
  
  comment_keys = openascii_comment_keys;
  keys = []; alt_name =[];
  
  /*   Read the file for cheking numberof of comment line
       if there is format definition lines
       the file size  will be defined on the createad function  */
  size_comment = strlen(comment_string);
  size_format_def = strlen(format_def_string);
  B_start = bookmark(f);
  
  for (tabpos=1; (line=rdline(f)) && line=="";tabpos++) ;
  if (line==string()) {
    if (prompt) write, "The file seems to be void or end of file";
    return [];
  }
  backup, f;

  if (!is_void(tab) && tab>1) {
    itab=1;
    do {
      for ( ; (line=rdline(f)) && line!="";tabpos++) ;
    if (line==string()) error, "I didn't find a table "+pr1(tab);
    backup, f;
    for ( ; (line=rdline(f)) && line==""; tabpos++) ;
    if (line==string()) error, "I didn't find a table "+pr1(tab);
    backup, f;
    itab++;
    } while(itab<tab);
    
  }
  tabpos--;
  
  size_c = 0; header=[];
  // Search number of comment lines
  do {
    linef=strtrim(rdline(f));
    cont=0;
    if (linef) {
      
      if (strpart(linef, 1:size_comment)==comment_string) {
        header=linef;
        size_c++;cont=1;
      }
      else if (strpart(linef, 1:size_format_def)==format_def_string) {
        if (!onlystruct && !noheader)
        openascii_read_format_line, strpart(linef,size_format_def+1:0) , keys, alt_name;
        size_c++;cont=1; 
      }
    }
  } while(cont);
  //for (size=1; strtrim(rdline(f)); size++);
  //close, f;
  backup, f, B_start ; backup,f;
  if (header)  header=strpart(strtrim(header),size_comment+1:0);
  else noheader = 1;
  no_op = no_op || noheader;
  if (noheader && !is_void(structure)) onlystruct=1;
  
  if (is_void(file_format)) file_format=_default_ASCII_format_file;

  
  Nadd = numberof(add_kw);
  for(i=1;i<=Nadd;i++) openascii_read_format_line, add_kw(i),keys,alt_name;
  
  // Read the file with format definition
  if (!onlystruct && !noheader && !is_void(file_format) && strtrim(file_format)) {
    g = open(file_format);
    while (line=rdline(g)) {
      openascii_read_format_line, line, keys, alt_name;
    }
    close, g;
  }
  else if (onlystruct) {
    openascii_get_struct_members, structure, name, type;
    
    for (i=1;i<=numberof(name);i++)
      openascii_read_format_line, name(i)+" "+type(i), keys, alt_name;
    
    if (noheader) {
      header = openascii_strjoin(name, " ");
      if (prompt) {tmptext = openascii_strjoin("-"(-:1:strlen(header)));
        write, format="%s\n\n","*** No header use structure definition to create it ***\n"+tmptext+"\n"+header+"\n"+tmptext;
      } 
      noheader=0;
    }
  }
  if (noheader && !is_void(strgrep)) {
    /* If there is no header Open ascii try to create one
     with the first data line */
    tabline = array(string, _max_key_in_header);
    sread, linef, tabline;
    tabline = tabline(where(tabline));
    tsize = numberof(tabline);
    if (is_void(valtype)) {
      test = strgrep("[a-zA-Z/;'\"?><,=%$#@!~`*:)(|\\}{]", tabline)(2,)>-1;
      type = merge2("string", "double",test);
    }
    else {
      if (typeof(valtype)=="struct_definition") {
        valtype=nameof(valtype);
      }
      else if(typeof(valtype)!="string") error, "valtype must be a string or a structure definition";
      type=valtype(-:1:ysize);
    }
    name = swrite(indgen(tsize),format="X%d");
    header = openascii_strjoin(name, " ");
    keys = array(KEYS,tsize);
    alt_name = array(ALT_NAME, tsize);
    keys.name=name; keys.type=type;  keys.num=indgen(tsize);
    alt_name.pos =1 ; alt_name.num = indgen(tsize) ; alt_name.name =name;
    if (is_void(valtype)) { 
      if (prompt) write, format="%s\n\n", "**** No header and no structure definition, I try to create one\n     with the first data line : all numerical values are \"double\" ***";
    } else { 
        if (prompt)  write, format="%s\n\n", "**** No header and no structure definition, you have specified a type \""+valtype+"\" I try to read the file with this type ****";
    }
  }
 
  
  if (!numberof(keys)) error, "I didn't find definitions to read this file";
  /* Check if there is special keys */
  if (numberof( (ou=where(keys(alt_name.num).type=="special")) )){
    keys(alt_name(ou).num).type = "string";
    _, comment_keys, alt_name(ou).name;
  }
    
  
  //extern _f_openASCII
    //_f_openASCII = open(file);
  // Go to the last comment line;
  //for (i=1; i<=size_c; i++) {header=rdline(_f_openASCII);}
  
  if (sep) {
      while( (p=strfind(sep,header))(0)!=-1){
      header = streplace(header, p, " ");
    }
  }
  
  header_tab = array(string, _max_key_in_header);
  sread, header, header_tab;
  header_tab = header_tab(where(header_tab));

  size_header = numberof(header_tab);
 
  header = array(ALT_NAME, size_header);
  header.name = header_tab;
  header.pos = indgen(size_header);
  operation = [];
  treated = array(int,size_header);

  // Now Do some operation
  if (!no_op) {
  do {
    cont=0;
    for (i=1;i<=size_header;i++) {
      if (!treated(i) && (op=vstrmatch(header(i).name, ["-","+","*","/"])) ) {

        tmptok = strtok(header(i).name, op);
        arg1 = tmptok(1);
        arg2 = tmptok(2);
        
        ou1 = where(arg1==header.name);
        ou2 = where(arg2==header.name);

        Nop = numberof(operation);
        /* Do not do operation if the two variables are in the header */
        if ( !(numberof(ou1) && numberof(ou2)) ) {
          
          if (numberof(ou1) && !numberof(ou2)) {
          grow, operation, openascii_treat_op(1,op,header(ou1).pos,header(i).pos,arg2,
                                    Nop+1);
          treated(i)=1; cont=1; 
        }
        else if (numberof(ou2) && !numberof(ou1)) {
          grow, operation, openascii_treat_op(2,op,header(ou2).pos,header(i).pos,arg1,
                                   Nop+1);
          treated(i)=1; cont=1; 
        }
        //else if (!numberof(ou1) && !numberof(ou2)) {}
        else if (Nop) {//If operation 
          ou1 = where(arg1==operation.name);
          ou2 = where(arg2==operation.name);
          if ( !(numberof(ou1) && numberof(ou2)) ) {
          if (numberof(ou1) && !numberof(ou2)) {
            grow, operation, openascii_treat_op(1,op,-operation(ou1).pos,header(i).pos,arg2, Nop+1);
            cont=1; 
          }
          if (numberof(ou2) && !numberof(ou1)) {
            grow, operation, openascii_treat_op(2,op,-operation(ou2).pos,header(i).pos,arg1, Nop+1);
            cont=1; 
          }
          }
        }
        }
      }
   
      
    } 
  } while(cont);
  }
  ou_good = ou_good_op =[];
  //check where alt_name = header;
  ou = openascii_comp_vecto (alt_name.name, header.name , test_alt, test_head);

  //check if there is same keys name  and take the max position in header
  oudb = openascii_op_doublon(keys(alt_name.num(ou(,1))).name,
                     -header(ou(,2)).pos, -header(ou(,2)).pos);   
   ou_good = ou(oudb,);
  
  //ou_good = ou;
  // Truncate header in fonction;
  bad_header = header(where(!test_head));
    
  if (numberof(operation)) {
    ou = openascii_comp_vecto (alt_name.name, operation.name , test_alt, test_head);
    if (numberof(ou)) {
      //oudb_op=openascii_op_doublon(alt_name.name(ou(,1)), alt_name.pos(ou(,1)), ,mnx);  
      //ou_good_op = ou(oudb,);
    }
    ou_good_op = ou;
  }
  
  //text_eval = "/* This file was created by OpenASCII the "+
  //  getdate()+" at "+gettime()+"*/\n";
  text_eval = "/* This file was created by OpenASCII */\n";
  //structure definition
  if (is_void(structure)) {
    /* Define a tmp structure  */
    if (is_void(struct_name))
      def_struct_name = "OPEN_ASCII_TMP";
    else
      def_struct_name = struct_name;
    
    text_eval += "\nstruct "+def_struct_name+" {\n"+
      openascii_strjoin(swrite (keys(alt_name(ou_good(,1)).num).type,
                       keys(alt_name(ou_good(,1)).num).name,
                       format="  %s  %s;\n"));
    if (numberof(ou_good_op)) {
      text_eval +=  openascii_strjoin(swrite(keys(alt_name(ou_good_op(,1)).num).type,
                                    keys(alt_name(ou_good_op(,1)).num).name,
                                    format="  %s  %s;\n"));
    }
    text_eval +="  }\n";
  
  
  }
  else { /* If user post a structure. */
    openascii_get_struct_members, structure , st_name, st_type;
    ou = openascii_comp_vecto (keys(alt_name(ou_good(,1)).num).name,
                      st_name, test_keys, test_st);
    if (!numberof(ou)) error, "any structure extension correspond to a header keys";
    ou_good = ou_good(ou(,1),);
    test_comp = openascii_check_compatibility(keys(alt_name(ou_good(,1)).num).type,
                                     st_type(ou(,2)));
    if (anyof(test_comp))
      ou_good = ou_good(where(test_comp), );
    else ou_good=[];
    
    if (numberof(ou_good_op)) {
      ou = openascii_comp_vecto (keys(alt_name(ou_good_op(,1)).num).name,
                      st_name, test_keys, test_st);
      if (numberof(ou)) {
       ou_good_op = ou_good_op(ou(,1) ,);
       test_comp = openascii_check_compatibility(keys(alt_name(ou_good_op(,1)).num).type,st_type(ou(,2)));
       if (anyof(test_comp))
         ou_good_op = ou_good_op(where(test_comp), );
       else ou_good_op=[];
      }
      else ou_good_op = [];
    }
    if (is_void(ou_good_op)&&is_void(ou_good))
      error, "There are No valid keys in your structure definition";
    def_struct_name = nameof(structure);
  }


  /* Variable declaration */
  sup = "var_";
  var_names = swrite( header.pos,  format=sup+"%d");
  var_type  = array("string", numberof(header));
  oudb = openascii_op_doublon(alt_name.name(ou_good(,1)), alt_name.pos(ou_good(,1)),
                     alt_name.num(ou_good(,1)));   
  var_type (ou_good(oudb,2)) = keys(alt_name.num(ou_good(oudb,1))).type;

  if (numberof(operation)) {
    ou = where(operation.col1>0);
    if (numberof(ou))
      var_type (operation.col1(ou)) = "double";
    ou = where(operation.col2>0);
    if (numberof(ou))
      var_type (operation.col2(ou)) = "double";
    
  }
  
  //define the function
  text_eval += "func "+func_name+" (&f,N, prompt=) {\n";
  //text_eval += swrite(file,format = "  file = \"%s\";\n");
  text_eval += "  if (is_void(prompt)) prompt=_PROMPT;\n";
  //text_eval += "  f = open(file);\n";
  text_eval += "  B_head = bookmark(f);\n";
  if (tabpos) {
    if (!is_void(tab))
      text_eval += openascii_put_write("I go to the table "+pr1(tab));
    text_eval += swrite(format="  for (i=1;i<=%d;i++) rdline,f;\n",tabpos);
  }
  text_eval += openascii_put_write("I Check number of lines to read ... ",,1);
  
  text_eval += swrite( size_comment,comment_string,
                      size_format_def,format_def_string, 
                      format="  cont=1;\n  for (Nc=-1; cont && (line=strtrim(rdline(f)));Nc++) {\n    if ( strpart(line,1:%d)==\"%s\" || strpart(line,1:%d)==\"%s\") cont=1;\n    else cont=0;\n }\n");
  text_eval += "  B_data = bookmark(f);\n";
  text_eval += "  if (is_void(N))";
  text_eval += " for (N=1; (line=rdline(f)) && line!=\"\"; N++) /* check file size */;\n";
  //text_eval += "  close, f;";
  text_eval += "  if (prompt) write ,Nc, N, format=\"%d comments + %d data lines\\n\";\n";
  
  text_eval += openascii_strjoin(swrite( var_names, var_type,
                                header.name, 
                                format="  %s =array(%s,N);// %s"), "\n")+"\n";
  
  if (sep)
     text_eval += swrite( "var_sep", "string",
                                " to read separations", 
                          format="  %s =array(%s,N);// %s")+"\n";
  // Sort with column number
  header = header(sort(header.pos));
  order = msort(header(ou_good(,2)).pos);
  ou_good = ou_good(order, );
  test_comment =  !is_void(comment_keys) &&
    anyof( header(0).name == comment_keys);
  
  //Read the file;
  
  text_eval +=  "\n  backup, f, B_data;\n  backup, f;\n";
  text_eval+=openascii_put_write("I Read the file ", "print(f)(2)+\"...\"",1);
  //text_eval +=  "\n  f = open(file);\n";
  //if (tabpos) {
  //  text_eval += "  /* go to the table "+pr1((tab?tab:1))+" */\n";
  //  text_eval += swrite(format="  for (i=1;i<=%d;i++) rdline,f;\n",tabpos);
  //}
  //text_eval += "  /* Skip comment lines */\n";
  //text_eval +="  for (i=1; i<=Nc; i++) { rdline,f;}\n";
  
  if (test_comment) {
    text_eval += "  /*Because last key is comment*/\n  var_last = array(string,N);\n";
    if (sep) {
      val_to_read =
        openascii_strjoin(
                          var_names(1:-2)+
                          merge2(",var_sep","",var_type(1:-2)!="string") ,
                          ",")+", var_last";
      format      = openascii_strjoin(openascii_type_to_format(var_type(1:-1), sep)(1:-1))+"%[^\\n]";
    }
    else {
      val_to_read = openascii_strjoin( var_names(1:-2), ",")+", var_last";
    format      = openascii_strjoin(openascii_type_to_format(var_type(1:-2)))+"%[^\\n]";
    }
  }
  else {
    if (sep)
      val_to_read =
        openascii_strjoin( var_names+merge2(",var_sep","",var_type!="string"&numberof(var_type)!=indgen(numberof(var_type))), ",");
    else 
      val_to_read = openascii_strjoin( var_names, ",");
    format  = openascii_strjoin(openascii_type_to_format(var_type, sep));
  }
  text_eval += "  statu =read( f, "+val_to_read+", format=\""+format+"\");\n\n";

  text_eval += "  /* Check if there is other tabs */\n";
  text_eval += "  while ((line=rdline(f)) && line==\"\") {}\n";
  text_eval += "  check_tab = line!=string();\n";
  
  //Close the file
  //text_eval += "  close , f;\n";
  if (test_comment) {
    text_eval += "  /* Because last key is comment */\n";
    text_eval += "  sread, var_last+\" _\", "+var_names(-1)+","+ var_names(0)+
      ",format=\""+ openascii_type_to_format(var_type(-1))+"%[^\\n]\";\n";
    Nvars = numberof(var_names)-1;
    text_eval += swrite(var_names(0),var_names(0), format="  %s = strpart(%s, 1:-2);\n");
    text_eval += "  /* --------------------------- */\n";
  }
  else 
    Nvars = numberof(var_names);
  if (sep) Nvars= numberof(var_names)+numberof(where(var_type!="string"))-
    (var_type(0)=="string" ? 0 : 1)-(test_comment ? 2 : 0);
  text_eval += "\n\n";
  text_eval += swrite( Nvars, format="  diff_statu = N*%d-statu;\n");
  warning_line = "  if ((diff_statu>0) & prompt) write, diff_statu, format=\"\\n WARNING %d elemements were not read !!!!\\n Check format types, N values or if there are obscure lines at the end of table!!!\\n\";\n";
  text_eval += warning_line;
  text_eval += "  else if (prompt) write, format=\"%s\\n\\n\", \"OK\";\n";
  //Operation :
  if (numberof(operation)) {
    text_eval += "  /* OPEARTION */\n";
    text_eval+=openascii_put_write("OPERATION -- OPERATION -- OPERATION --");
    op_names = swrite(indgen(numberof(operation)), format="ope_%d");
    for (i=1;i<=numberof(operation);i++) {
      opname1 = (operation(i).col1>0 ? header(operation(i).col1).name :
                 operation(-operation(i).col1).name);
      opvar1 = (operation(i).col1>0 ?swrite(operation(i).col1, format=sup+"%d") :
                op_names((-operation(i).col1)) );

      opname2 = (operation(i).col2>0 ? header(operation(i).col2).name :
                 operation(-operation(i).col2).name);
      opvar2 = (operation(i).col2>0 ?swrite(operation(i).col2, format=sup+"%d") :
                op_names(-operation(i).col2) );
      
      
      text_eval += openascii_put_write(swrite(opname1, 
                                              operation(i).op,
                                              opname2, 
                                              operation(i).name,
                                              format="(%s) %s (%s) => %s"));
      text_eval += swrite(op_names(i), opvar1,
                          operation(i).op, opvar2,
                          format="  %s = %s %s %s ;\n\n");
    }
  }

  text_eval += "\n/*   Fill the structure */ \n";
  text_eval += "  out = array("+def_struct_name+",N);\n\n";
  
  // Fill the structure.
  test_c = openascii_check_compatibility( keys(alt_name(ou_good(,1)).num).type,
                                var_type(ou_good(,2)) );
  ou_c = where(test_c);
  textc1 = merge2("", "/*",test_c);
  textc2 = merge2("", "Invalid types */",test_c);
  text_eval += openascii_strjoin(swrite(textc1,
                               keys(alt_name(ou_good(,1)).num).name,
                               var_names(ou_good(,2)),textc2,
                              format=" %s out.%s = %s; %s\n"))+"\n";
  if (numberof(ou_good_op)) {
    test_c = openascii_check_compatibility( keys(alt_name(ou_good_op(,1)).num).type,
                                  "double" );
    textc1 = merge2("", "/*",test_c);
    textc2 = merge2("", "Invalid types */",test_c);
    text_eval += openascii_strjoin(swrite(textc1,
                                 keys(alt_name(ou_good_op(,1)).num).name,
                                 op_names(ou_good_op(,2)),
                                 textc2,
                                 format=" %s out.%s = %s; %s\n"))+"\n";
    
  } 
  // free som memories
  //text_eval += openascii_strjoin(var_names, "=")+"= [];\n";
  //if (numberof(operation))
  //   text_eval += openascii_strjoin(op_names, "=")+"= [];\n";


  
  
  if (numberof(bad_header)) {
    text_eval += openascii_put_write ("\\nBAD KEY -- BAD KEY -- BAD KEY")
      text_eval += openascii_put_write (openascii_strjoin(swrite( bad_header.name, format="%s  "))+"\\n");
  }
    
  text_eval += openascii_put_write("\\nGOOD KEY -- GOOD KEY -- GOOD KEY");
  text_eval += openascii_put_write(openascii_strjoin(swrite (header(ou_good(,2)).name , keys(alt_name(ou_good(,1)).num).name, strpart(keys(alt_name(ou_good(,1)).num).type,1:1), format="%s => %s (%s)"), " ** "));
  
  if (numberof(ou_good_op))
    text_eval += openascii_put_write(openascii_strjoin(swrite ( operation(ou_good_op(,2)).name , keys(alt_name(ou_good_op(,1)).num).name, format="%s => %s"), " ** ")+"\\n");    


  
  if (warning_line)
    text_eval += warning_line;

  text_eval += "  if (check_tab) {\n";
  text_eval += "  backup, f;\n";
  text_eval +=  openascii_put_write("\\nThere is a following tab in this file.");
  text_eval += " }\n"// else { close, f;}";
  //cloe the function.
  text_eval += "\n  return out;\n}\n";
  text_eval += "tmp_func_openASCII = "+func_name;
  openascii_eval, text_eval, tmp=code_file , debug=debug;
  
  st_out = tmp_func_openASCII(f,N,prompt=prompt);
  
  return st_out;
}
OpenASCII=OpenASCII2;
ReadASCII=ReadASCII2;

func openascii_treat_op (Narg , op, num1, num2, name, pos) {
/* DOCUMENT Fonction used by OpenASCII
   SEE ALSO:
 */
  if (op=="+")
    return OPERATION(op="-", col1=num2, col2=num1,name=name,pos=pos);
  else if (op=="-" && Narg==1) 
    return OPERATION(op="-", col1=num1,col2=num2,name=name,pos=pos);
  else if (op=="-" && Narg==2) 
    return OPERATION(op="+", col1=num1,col2=num2,name=name,pos=pos);
  else if (op=="/" && Narg==1)
    return OPERATION(op="/", col1=num1,col2=num2,name=name,pos=pos);
  else if (op=="/" && Narg==2)
    return OPERATION(op="*", col1=num1,col2=num2,name=name,pos=pos);
  else if (op=="*")
    return OPERATION(op="/", col1=num2, col2=num1,name=name,pos=pos); 
}

func vstrmatch (s,v) {
  size = numberof(v);
  for (i=1; i<=size; i++) if (strmatch(s,v(i))) return v(i);
  return 0;
}


func openascii_op_doublon (value, rule, rule2, eps=) {
/* DOCUMENT Function used by OpenASCII
   return index of doublon. Doublon are choosen function tu rule and
   rule2
   SEE ALSO:
 */
  if (is_void(eps)) eps = 1e-20;
  N = indgen(numberof(value));
  num=[];
  structure = structof(value(1));
  do {
    if (structure=="double")
      test = abs(value-value(1))<eps;
    else
      test = value==value(1);
    ou = where(test);

    ou = ou(msort(rule(ou),rule2(ou)));
    //grow, num, N(ou(rule(ou)(op)));
    grow, num, N(ou(1));
    N = N(where(!test));
    value =   value(where(!test));
    rule =    rule(where(!test));
    rule2 =   rule2(where(!test));
  } while(!is_void(value));
  return num;
}


func openascii_get_struct_members (st, &name, &type) {
/* DOCUMENT openascii_get_struct_members (st, &name, &type)
   Function used by OpenASCII.
   return names and types of the strucure definition
   SEE ALSO:
 */
  if (typeof(st)!="struct_definition")
    error,  "expexting structure definition";
  st_text = strtrim(print(st)(2:-1));
  tab = strtok(st_text, " ");
  name= strpart(tab(2,) , 1:-1); type=tab(1,);
}

func openascii_check_compatibility (var1, var2) {
/* DOCUMENT Function used by OpenASCII
   check if var1 and var2 have a compatible type;  
   SEE ALSO:
 */
  var1 = merge2 ("numerical", var1, (var1=="int"| var1=="long"|var1=="double"|
                              var1=="float"|var1=="char"));
  var2 = merge2 ("numerical", var2, (var2=="int"| var2=="long"|var2=="double"|
                              var2=="float"|var2=="char"));
  return var1==var2;
}

func openascii_put_write (text,val, nopass) {
/* DOCUMENT Function used by OpenASCII
   SEE ALSO:
 */
  return "  if (prompt) write, format=\"%s"+(nopass ? "" : "\\n")+
    "\", \""+text+"\""+(val ? "+"+val : "")+";\n";
}
func openascii_strjoin(str, glue)
/* DOCUMENT openascii_strjoin(str)
       -or- openascii_strjoin(str, glue)
     Join strings from array STR into a single long string.  The string GLUE
     (default "") is used between each pair of element from STR.

   SEE ALSO strcut */
{
  if ((n= numberof(str)) >= 1) {
    s= str(1);
    if (glue) for (i=2 ; i<=n ; ++i) s+= glue+str(i);
    else      for (i=2 ; i<=n ; ++i) s+= str(i);
    return s;
  }
  return "";
}

func openascii_comp_vecto (X1, X2, &test1 , &test2, determinist=, same=)
/* DOCUMENT Function used by OpenASCII
   Check pairs in to tabular
   SEE ALSO:
 */
{
  size1 = numberof(X1); size2 = numberof(X2);
  
  order1 = (size1>1 ? msort(X1) : [1]);
  if(same) order2=order1;
  else order2 = (size2>1 ? msort(X2) : [1]);
  
  X1 = X1(order1);  X2 = X2(order2);

  test1 = array(int, numberof(X1));
  test2 = array(int, numberof(X2));
  
  
  match1 = match2 = [];
  for (i = j = 1; (i <= size1) && (j <= size2);) {
    if (X1(i) < X2(j)) i++;
    else if (X1(i) > X2(j)) j++;
    else {
      test = 1;
      first_j = j;
      for (j = first_j; ( test ) && (same ? j<i : (j <= size2)); j++) {
        if ( (test = (X1(i) == X2(j)) )) {
          _, match1, order1(i);
          _, match2, order2(j);
          test1(order1(i)) = 1;
          test2(order2(j)) = 1;
        }
      }
      j = first_j;
      i++;
      //i = (i==size1 ? i : i++);
    }
  }
  if (!numberof(match1)) {
    return [];
  }
  return [match1, match2];
}

func openascii_type_to_format (type, sep) {
/* DOCUMENT Function used by OpenASCII
   Change type to format for read or write function
   SEE ALSO:
 */
  size = numberof(type);
  def_type = ["long", "int", "double", "float", "string" , "char"];
  def_format = ["%d", "%d" , "%f"   ,  "%f"   , "%s"   ,   "%d"];

  test = array(int, dimsof(type));
  type_size = numberof(def_type);
  for (i=1;i<=type_size;i++)
    test = test+(type==def_type(i))*i;
  if (nallof(test)) error, "There are invalid type";
  out_format = def_format(test);
  if (sep && numberof(out_format)>1) {
    sep = "%[^"+strpart(sep, 1:1)+"]"+strpart(sep, 1:1);
    out_format(1:-1) =  merge2(out_format(1:-1)+sep, sep, out_format(1:-1)!="%s");
  }
  return out_format;
}



func openascii_read_format_line (line, &keys, &alt_name) {
/* DOCUMENT Function used by OpenASCII
   Read a line definition of OpenASCII and fill keys and alt_name vars 
   SEE ALSO:
 */
  size = numberof(line);

  if(!is_void(keys)) num=max(keys.num)+1;
  else num=1;
  tmp_v = array(string,_max_alt_key_number);
  sread, line, tmp_v;    
  tmp_v = tmp_v(where(tmp_v));
  name = tmp_v(1);
  if (numberof(keys)) {
    if (anyof(name==keys.name))
      return ;
  }
  grow, keys , KEYS(num=num, 
                    name=name,
                    type=tmp_v(2));
  if (numberof(tmp_v)<=2) {
    grow, alt_name,ALT_NAME(name=name,num=num,pos=1);
    return ;
  }
  tmp_size = numberof(tmp_v(3:0));
  tmp_alt_name = array(ALT_NAME,tmp_size);
  tmp_alt_name.name = tmp_v(3:0);
  tmp_alt_name.num  = num;
  tmp_alt_name.pos  = indgen(tmp_size);
  grow, alt_name, tmp_alt_name;
  
}
  

func openascii_eval(code, tmp=, debug=)
/* DOCUMENT openascii_eval, code;
       -or- openascii_eval (code);
       
     This is a copy of eval by Eric Thibau (utils.i)
     Include here for OpenASCII.
     
     Evaluate CODE given as a string  or as an array of strings (considered
     as  different lines  in the  script).  Since  CODE can  be dynamically
     build,   this  routine   allows  the   execution  of   virtually  (see
     hints/restrictions below)  any Yorick's code  (e.g. dynamic definition
     of  structures,  of functions,  etc.).   For  instance, the  following
     statement defines a new structure:
       eval, "struct NewStruct {string a; long b; float c, d;}";

     Since  the script  gets evaluated  at the  scope level  of  the "eval"
     routine some local variables of the  "eval" routine may be used in the
     script:
       "eval_tmp"    contains  the  name of  the temporary script  file and
                     must not be changed by the script;
       "eval_debug"  contains the value of  the keyword DEBUG and  must not
                     be changed by the script;
       "eval_code"   contains the value of the argument CODE;
       "eval_result" is returned by "eval", its contents may be defined into
                     the script.
     Note: impredictible  results may  occur if CODE  changes the  value of
     symbols "eval_tmp" and "eval_debug".

     Keyword TMP  can be  used to  specify the file  name of  the temporary
     script.  The default file name is:
       "$YORICK_EVAL_TMP"      if environment variable "YORICK_EVAL_TMP" set;
                               is set;
       "/tmp/$USER-eval_tmp.i" if environment variable "USER" set;
       "~/.eval_tmp.i"         otherwise.

     If  keyword DEBUG  is true  (non-zero and  non-nil), the  name  of the
     temporary file is printed out and the file is not removed.


   SEE ALSO: include. */
{
  /* The "eval_" prefix is used in order to somewhat protect local variables
     from caller's code. */
  local eval_result, eval_tmp, eval_code, eval_debug;
  eq_nocopy, eval_tmp,   tmp;
  eq_nocopy, eval_debug, debug;
  eq_nocopy, eval_code,  code;

  /* Dump script into a temporary file, then include it. */
  if (is_void(eval_tmp)) {
    /* Create default name for yorick temporary code. */
    eval_tmp = get_env("YORICK_EVAL_TMP");
    if (! eval_tmp) {
      eval_tmp = get_env("USER");
      eval_tmp = (eval_tmp ? "/tmp/"+eval_tmp+"-" : "~/.") + "eval_tmp.i";
    }
  }
  write, format="%s\n", open(eval_tmp, "w"), eval_code;
  include, eval_tmp, 1;
  if (eval_debug) {
    write, format="Yorick code written in \"%s\"\n", eval_tmp;
  } else {
    remove, eval_tmp;
  }
  return eval_result;
}
