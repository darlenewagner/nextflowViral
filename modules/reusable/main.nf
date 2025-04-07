process checkExecutables {

  input:
  val executable
  
  output:
  stdout

  script:
"""
   if type "${executable}" &> /dev/null; then
      echo "true"
  else
      echo "false"
  fi  
"""    
}


