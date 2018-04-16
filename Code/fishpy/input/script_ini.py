import os 

for sensi in [0,1,2]:
    for fsky in [0.1,0.2,0.4]:
        for lensing in [True,False]:
            for rayleigh in [True,False]:
                with open('SO.ini','r') as f:
                    newtext = ''
                    update_f = True
                    update_ray = True
                    for line in f:
                        if 'sensitivity' in line : 
                            line1,line2 = line.split('=')
                            line2 = ' = {:d} \n'.format(sensi)
                            line = line1 + line2
                        if 'fsky' in line : 
                            if update_f:
                                line1,line2 = line.split('=')
                                line2 = ' = {:3.1f} \n'.format(fsky)
                                line = line1 + line2
                                update_f = False
                            else:
                                pass
                        if 'lensing' in line:
                            line1,line2 = line.split('=')
                            line2 = ' = {} \n'.format(lensing)
                            line = line1 + line2
                        if 'rayleigh' in line:
                            if update_ray:
                                line1,line2 = line.split('=')
                                line2 = ' = {} \n'.format(rayleigh)
                                line = line1 + line2
                                update_ray = False
                            else:
                                pass
                        newtext += line
                str1 = 'l' if lensing else 'nol'
                str2 = 'r' if rayleigh else 'nor'
                filename = 'SO_{:d}_{:3.1f}_{}_{}.ini'.format(sensi,fsky,str1,str2)
                with open('SO_ini/'+filename,'w') as f:
                    f.write(newtext)
                print('Done for {:d} {:3.1f} {} {}'.format(sensi,fsky,str1,str2))
                        
