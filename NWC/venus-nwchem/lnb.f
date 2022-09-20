c
c   This function returns the location of last non-blank character.
c   The usage is chr(1:lnb(chr)) will give 'some' instead of 'some   '
c   This routine works fine with linux, but have some problem on SGI
c
        integer function lnb(ch)
        character*(*) ch
        do i = len(ch), 1, -1
            if (ch(i:i) .ne. ' ') then
                 lnb = i
                 return
            endif
        enddo
        lnb = 0   ! error ???
        end
