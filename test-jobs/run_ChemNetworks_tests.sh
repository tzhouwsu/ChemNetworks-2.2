#!/bin/bash

echo "here are test jobs for ChemNetworks-2.2"
echo ""

echo "test1.."
cd ./test1/
cp ../../ChemNetworks-2.2.exe ./
./ChemNetworks-2.2.exe Input-test1 water1.xyz
lab=$?

if test $lab -eq 0; then
  echo "test1 passed"
  echo ""

  echo "test2.."
  cd ../test2/
  cp ../../ChemNetworks-2.2.exe ./
  ./ChemNetworks-2.2.exe Input-test2 water1.xyz solB1.xyz
  lab=$?
  if test $lab -eq 0; then
     echo "test2 passed"
     echo ""

     echo "test3.."
     cd ../test3/
     cp ../../ChemNetworks-2.2.exe ./
     ./ChemNetworks-2.2.exe Input-test3 water1.xyz solB1.xyz solC1.xyz
     lab=$?
     if test $lab -eq 0; then
        echo "test3 passed"
        echo ""

        echo "test4.."
        cd ../test4/
        cp ../../ChemNetworks-2.2.exe ./
        ./ChemNetworks-2.2.exe Input-test4 water1.xyz solC1.xyz solB1.xyz
        lab=$?
        if test $lab -eq 0; then
           echo "test4 passed"
           echo ""
        else
           echo "test4: Error"  # if test4 failed
        fi

     else               # if test3 failed
        echo "test3: Error!"
     fi 

  else               # if test2 failed
     echo "test2: Error!"
  fi 

else       # if test1 failed
  echo "test1: Error!"
fi



