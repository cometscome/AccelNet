from pathlib import Path
import shutil
root=Path('/tmp/accelnet-angular-contraction')
for name in ['source','source-pairs']:
 dst=root/(name+'-local');shutil.copytree(root/name,dst,dirs_exist_ok=True)
 p=dst/'AccelNetPredictor/src/accelnet_target_kernels.f90';s=p.read_text();a=s.index('    subroutine generic_values');end=s.index('    subroutine generic_contract',a);v=s[a:end]
 v=v.replace('kind,t1,t2,group,bb','kind,t1,t2,group,bb,degree,d')
 v=v.replace('cosine,total12,v12,dv12\n','cosine,total12,v12,dv12,t,power,hvalues(0:16)\n')
 v=v.replace('cosine,total12,v12,dv12)', 'cosine,total12,v12,dv12,t,power,hvalues)')
 b=v.index('                if (features(10,b,s) == b) then');e=v.index('                group =',b)
 v=v[:b]+'''                degree = features(12,b,s)
                hvalues = 0
'''+v[e:]
 b=v.index('                            bb = b');e=v.index('                        end do',b)
 v=v[:b]+'''                            if (degree > 0) then
                                t = 0.5_real64*(1+fp(4,b,s)*cosine); power = 1
                                do d=1,degree
                                    power = power*t
                                    hvalues(d) = hvalues(d)+v*power
                                end do
                            else
                                total = total+v*angular_value(cosine,fp(4,b,s),fp(5,b,s),features(5,b,s))
                            end if
'''+v[e:]
 v=v.replace('                if (features(10,b,s) == 0) g(row,b) = total','''                if (features(10,b,s) == b) then
                    bb = b
                    do while (bb /= 0)
                        if (degree > 0) then
                            g(row,bb) = hvalues(features(5,bb,s))
                        else
                            g(row,bb) = total
                        end if
                        bb = features(11,bb,s)
                    end do
                else
                    g(row,b) = total
                end if''')
 s=s[:a]+v+s[end:]
 lines=[]
 for line in s.splitlines():
  while len(line)>130 and (not line.lstrip().startswith('!') or line.lstrip().startswith('!$omp')):
   k=line.rfind(',',0,110)
   if k<0:break
   lines.append(line[:k+1]+' &');line=('        !$omp& ' if line.lstrip().startswith('!$omp') else '            ')+line[k+1:].lstrip()
  lines.append(line)
 p.write_text('\n'.join(lines)+'\n')
