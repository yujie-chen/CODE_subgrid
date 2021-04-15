				implicit none

				integer :: MX, NY, Ratio,nt,it,I, &
                            beg_n,end_n,num_n, beg_e,end_e,num_e, &
                            beg_w,end_w,num_w, beg_s,end_s,num_s
				real       dt,time,ele,u,v,time_tot, &
                           A1,T1,w1,A2,T2,w2
				
				
				MX = 0
                                NY = 24
				Ratio = 1
				
				beg_n = 0
				num_n = 0
				end_n = 0

				beg_e = 0
                                num_e = 0
				end_e = 0
				
				beg_w = 1
                		num_w = NY
         		        end_w = NY

			        beg_s = 0
                                num_s = 0
                                end_s = 0

				



!                A1 = 0.75
                 A1 = 0.06
                 T1 = 355.0
                 w1 = 2*3.1415926/T1

               A2=0.0
               w2=0.0

               time_tot = T1*11
               dt = 10.0
               nt = ceiling(time_tot/dt)
                write(*,'(A10,I10)') 'nt = ',nt

!                A2 = 0.25
!                T2 = 12.42*3600
!                w2 = 2*3.1415926/T2
				
                open(1,file='bc_elevationa.txt')
				open(2,file='bc_tidea.txt')


				write(2,*)'boundary condition for Ideal vegetation root case'
				write(2,*)'boundary info: num of points, start point'				
				write(2,*)'EAST'
				write(2,*) num_e,beg_e
				write(2,*)'WEST'
				write(2,*) num_w,beg_w
				write(2,*)'SOUTH'
				write(2,*) num_s,beg_s
				write(2,*)'NORTH'
				write(2,*) num_n,beg_n
				write(2,*)'TIME SERIES'
				
				u = 0.0
				v = 0.0
				do it = 1, nt
					time = (it-1)*dt
 !                   ele = A1*sin(w1*time)+A2*sin(w2*time)
                     ele = A1*sin(w1*time)
                    !write(*,'(A10,I10)') 'it = ', it
                    !write(*,'(A10,F10.2)') 'time = ', time
                    !write(*,'(A10,F10.5)') 'ele = ', ele
                    write(1,'(F16.2,F16.5)') time, ele

					write(2,*)time
					
					write(2,*)'EAST SIDE'
 !                   write(2,100)(u,I=beg_e,num_e)
 !                   write(2,100)(u,I=beg_e,num_e)
 !                   write(2,100)(ele,I=beg_e,num_e)
					
					write(2,*)'WEST SIDE'
					
                    write(2,100)(u,I=beg_w,num_w)
                    write(2,100)(u,I=beg_w,num_w)
                    write(2,100)(ele,I=beg_w,num_w)

					write(2,*)'SOUTH SIDE'
					
  !                  write(2,100)(u,I=beg_s,num_s)
  !                  write(2,100)(u,I=beg_s,num_s)
  !                  write(2,100)(ele,I=beg_s,num_s)

					write(2,*)'NORTH SIDE'

   !                 write(2,100)(u,I=beg_n,num_n)
   !                 write(2,100)(u,I=beg_n,num_n)
   !                 write(2,100)(ele,I=beg_n,num_n)

				enddo
				
100     format(5E16.6)
																		
				end
