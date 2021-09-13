
#ifndef CALL_SPICE
#define CALL_SPICE

#include "config.h"
#include <sys/types.h>

#define LINE_LENGTH 1024

void pname(Configuration *);
pid_t start_spice(const Configuration *C, double accuracy, double *pc, double *po, const char *call,
                  const char *returnn);
void call_spice(const Configuration *C, double accuracy, double *pc, double *po, const char *call,
                const char *returnn);
int spice_dice(Configuration *);
int generic_spice_files(void);

/* generic spice file .malt.run */
#define MALT_RUN "\n\
*.malt.run: called by *.call\n\
*run simulation and print node_name vectors to a file\n\
.control\n\
\n\
source $pname\n\
if ($#param == 1)\n\
  $param\n\
end\n\
source $circuit\n\
\n\
run\n\
write $n_return $node_name\n\
\n\
set noaskquit\n\
\n\
.endc\n\
"

/* generic spice file .malt.binsearch */
/* ***** checking pc everytime is highly inefficent for corner vector margins analysis ***** */
/* the binary search has been modified to handle multiply/divide */
/* (instead of plus/minus) in log space */
#define MALT_BINSEARCH "\n\
*.malt.binsearch: called by *.call\n\
*find the operating boundary between the two points pc and po\n\
*if accuracy=0, see if pc parameter set is pass or fail\n\
.control\n\
\n\
*add codeblocks: these will have multiple accesses\n\
*can $circuit be a codeblock?\n\
codeblock $pname -a\n\
if ($#param == 1)\n\
  codeblock $param -a\n\
end\n\
if ($#passf == 1)\n\
  codeblock $passf -a\n\
end\n\
codeblock .malt.passfail -a\n\
\n\
source $envelope\n\
pegged=0\n\
\n\
*open the return file\n\
echo > $return\n\
\n\
*check inner point...except for corners\n\
param=pc\n\
if dashc == 1\n\
  failed=0\n\
else\n\
  .malt.passfail\n\
end\n\
echo $&failed >> $return\n\
\n\
* nominal vectors\n\
if dasht == 1\n\
  dasht = 2\n\
endif\n\
\n\
*if we failed, or if we are doing a point-check, bailout\n\
*po[0]=0 only for point-check\n\
\n\
*check outer point\n\
if (failed == 0 and po[0] <> 0)\n\
  param=po\n\
  .malt.passfail\n\
\n\
*find margin\n\
  if failed == 0\n\
    pegged=1\n\
  else\n\
    *this conditional-hack prevents division by zero for the case of pl == 0\n\
    *all because we are doing vector math\n\
    deltal=sqrt(((pl == 1)*po+(pl == 0))/((pl == 1)*pc+(pl == 0)))\n\
    delta =0.5*(po-pc)\n\
    param=(pl == 0)*(param-delta) + (pl == 1)*(param/deltal)\n\
    while (delta[0] > 1)\n\
      deltal=sqrt(deltal)\n\
      delta =0.5*delta\n\
      .malt.passfail\n\
      if failed=1\n\
        param=(pl == 0)*(param-delta) + (pl == 1)*(param/deltal)\n\
      else\n\
        param=(pl == 0)*(param+delta) + (pl == 1)*(param*deltal)\n\
      end\n\
    end\n\
  end\n\
  echo $&param >> $return\n\
\n\
  *trace marginal vectors\n\
  if dasht == 1 or dasht == 2\n\
    dasht=3\n\
    * margin point might pass or fail\n\
    .malt.passfail\n\
    if pegged == 0\n\
      if failed == 1\n\
        param=(pl == 0)*(param-delta) + (pl == 1)*(param/deltal)\n\
      else\n\
        param=(pl == 0)*(param+delta) + (pl == 1)*(param*deltal)\n\
      end\n\
    end\n\
    * this point ought to be the opposite (fail or pass)\n\
    if pegged == 1\n\
      pegged=2\n\
    end\n\
    .malt.passfail\n\
  end\n\
\n\
end\n\
set noaskquit\n\
quit\n\
\n\
.endc\n\
"

/* generic spice file .malt.passfail */
#define MALT_PASSFAIL "\n\
*.malt.passfail: called by .malt.binsearch\n\
.control\n\
\n\
*node names, envelope values, and step values contained in envelope file\n\
*the envelope file is sourced just once, in binsearch\n\
\n\
$pname\n\
if ($#param = 1)\n\
  $param\n\
end\n\
source $circuit\n\
failed=0\n\
\n\
*conserve memory\n\
save $node_name\n\
*step through time\n\
step_index = 0\n\
step_total = 0\n\
dowhile (step_index < length(step_value)) and failed=0\n\
  step $&step_value[$&step_index]\n\
  step_total_old = step_total\n\
  step_total     = step_total+$&step_value[$&step_index]\n\
\n\
  *echo step index: $&step_index\n\
\n\
  *step through list of nodes to check\n\
  node_index = 0\n\
  dowhile node_index < $#node_name\n\
    if ($node_name[$&node_index]) > ($node_hi[$&node_index])[0,length($node_name[$&node_index])-1]\n\
      failed=1\n\
      echo Node $node_name[$&node_index] failed hi envelope. Time step in range $&step_total_old - $&(step_total-1)\n\
    end\n\
    if ($node_name[$&node_index]) < ($node_lo[$&node_index])[0,length($node_name[$&node_index])-1]\n\
      failed=1\n\
      echo Node $node_name[$&node_index] failed lo envelope. Time step in range $&step_total_old - $&(step_total-1)\n\
    end\n\
    node_index=node_index+1\n\
  end\n\
  step_index=step_index+1\n\
end\n\
\n\
*execute manual passfail file\n\
if (failed = 0) and ($#passf = 1)\n\
  $passf\n\
end\n\
\n\
*trace marginal vectors\n\
if dasht=1 or dasht=3\n\
  if failed=1\n\
    * complete the simulation\n\
    while (step_index < length(step_value))\n\
      step $&step_value[$&step_index]\n\
      step_index=step_index+1\n\
    end\n\
  end\n\
endif\n\
if dasht=1\n\
  write $n_return $node_name\n\
endif\n\
if dasht=3\n\
  if failed=1 or pegged=2\n\
    write $f_return $node_name\n\
  else\n\
    write $p_return $node_name\n\
  end\n\
end\n\
\n\
free $node_name yes\n\
\n\
.endc\n\
.end\n\
"

#endif
