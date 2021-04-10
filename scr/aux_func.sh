message()
{
echo -e "\n"$1"\n"
}

confirm()
{
echo; read -p "Press Enter to proceed or CTRL-C to abort..."; echo
}

process_message()
{
echo -ne "\n"$1 ... ; cmd "${2}"; echo -e "done.""\n"
}

abort()
{
REASON="$1"

local i=0

while caller $i; do ((i++)); done

message "$REASON Abort."

exit 1
}

check_dir()
{
if [ ! -d $1 ]; then abort "$FUNCNAME: $1 <- no such directory."; fi

if [ $2 ]; then

   if [ $2 != "e" ]; then abort "$FUNCNAME: the second argument is not understood."; fi

   LIST="echo `ls $1`"

   if [ ! -n "$($LIST)" ]; then abort "$FUNCNAME: $1 <- the directory is empty."; fi

fi
}

check_file()
{
if [ ! -f $1 ]; then abort "$FUNCNAME: $1 <- no such file."; fi
}

check_file_dir()
{
if [ ! -f $1 ] && [ ! -d $1 ]; then abort "$FUNCNAME: $1 <- no such file or directory."; fi
}

check_if_exists()
{
if [ -z $2 ]; then abort "$FUNCNAME: error, no second argument."; fi

if [ $2 == "f" ];  then check_file $1;     fi

if [ $2 == "d" ];  then check_dir $1;      fi

if [ $2 == "fd" ]; then check_file_dir $1; fi
}

change_dir()
{
check_dir $1; cmd "cd $1"; message "DIRECTORY: $1"
}

exit_if_error()
{
if [ "$?" -ne 0 ]; then

   if [ ! "$1" ]; then abort "Abnormal termination: the exit status of the last command is not zero."; fi

   if [ "$1" ];   then abort "$1"; fi

fi
}

exit_if_noerror()
{
if [ "$?" -eq 0 ]; then abort "$1"; fi
}

cmd()
{
$1; exit_if_error "$2"
}

bg_cmd()
{
$1 &

exit_if_error
}

insert()
{
if [ ! "$3" ];                                     then abort "$FUNCNAME: error, no third argument.";                    fi
if [ "$3" ] && [ "$3" != "r" ] && [ "$3" != "a" ]; then abort "$FUNCNAME: error, the third argument is not understood."; fi
if [ "$4" ] && [ "$4" != "e" ];                    then abort "$FUNCNAME: error, the forth argument is not understood."; fi

if [ "$4" ]; then ECHO='echo -e'; else ECHO='echo'; fi

if [ "$3" == "r" ]; then $ECHO $1 >  $2; fi
if [ "$3" == "a" ]; then $ECHO $1 >> $2; fi
}

xterm_cmd()
{
insert "$1" "$XTERM_CMD_FILE" "a"
}

create_dir()
{
cmd "mkdir -p $1"
}

rm_dir()
{
if [ -d $1 ]; then chmod u+rwx -R $1; process_message "Directory $1 exists. Removing" "rm -r $1"; fi
}

rm_file()
{
if [ -f $1 ]; then cmd "rm $1"; fi
}

cp_file()
{
FILE=$1

DEST=$2

OPTION=$3

CP_CMD="cp"

if [ $OPTION ] && [ "$OPTION" == "v" ]; then CP_CMD="cp -v"; fi

if [ -f $FILE ]; then cmd "$CP_CMD $FILE $DEST"; fi
}

copy()
{
#for file in $1; do cp -v $file $2; done
for file in $1; do ln -s $file -t $2; done
}

#catenate()
#{
#ls $1/* > $2

#SORTED_NAMES=$(sort -V $2)

#rm $2

#for NAME in $SORTED_NAMES; do cat $NAME >> $3; done
#}

catenate()
{
if [ -d $1 ] && [ "$(ls -A $1)" ]; then

    sorted=$(ls $1/* | sort -V)

    for file in $sorted; do cat $file >> $2; done

fi
}
