FasdUAS 1.101.10   ��   ��    k             l      �� ��   y
Script to pipe selected Apple Mail messages to Yojimbo notes

Author: Steve Samuels (samplerx@[Protected])
File: yojimbo23.txt
PipeMailtoYojimbo.scpt

Version:23
Date: 2007/04/16

This script is based on one by Jim Correia of Barebones, cited 
below. It adds the following information to to the "Comments" of an 
imported email note. If the comments Column is visible in Yojimbo's 
Item List Pane, then the information can be viewed without opening 
the note.

* Added to Comments
1) the date sent (yyyy/mm/dd)
2) the name or address of sender or recipient
3) the first 60 characters of the message itself it, after 
removing returns, tabs, and most extra spaces.

* Adds tags: "email" & either "in" or "out"

After running the script, a Comment will look something like:

"2007/04/08 <: NYTimes.com---Your subscription is up soon. Please 
rene" (an incoming message)
or
"2007/04/09 >: Sam---Hi, Sam! I can't believe that you won the 
lotte" (an outgoing message)

Notes:

1. Text will not import for some rich text messages with 
attachments. The date sent, sender & tags are added, but 
AppleScript's "content of message" function returns only variants of 
"??". If this happens, open the original message & cut & paste the 
content into the note created by the script.

2. Graphics are removed and replaced by question marks. If you wish 
to retain the graphics, import the mail via the "Save PDF to Yojimbo" 
option in Mail's Print window.

3. A message is classified as outgoing if the sender's address 
matches one of the user's addresses. If an account has been deleted, 
there is no address to match, and all messages in the account will be 
classified as incoming. To fix this behavior, create and inactivate a 
dummy account. Only the email address must be correct; the other 
settings are irrelevant.

4. You can reset the number of characters that are added to comments 
from the body of the message. The setting statement appears just 
after the "on run" command at the start of the script proper.

5. If you have selected messages in more than one Message Viewer, the 
script transfers only those from the Viewer that was opened first.

6. Changing the font of Yojimbo's item list may make the the list 
more readable. (I now use Lucinda-Grande Bold 11)

Acknowledgements:

The original form for this script was based on MailtoYojimbo, by Jim 
Correia of Bare Bones Software, posted in the Yojimbo Talk Forum: 
http://www.listsearch.com/Yojimbo/Message/index.lasso?2169
& at //www.hawkwings.net/2006/09/04/script-to-pipe-emails-into-yojimbo/

I got the idea for adding extended message information to the 
comments field after seeing the "Summary" field in SOHO Notes.

The error dialog was copied from MailToYojimbo-v1.1 by Drummond 
Field, posted at
http://www.hawkwings.net/2007/02/05/email-to-yojimbo-script-with-pdf- 
support

The source for the minimum function was:
http://www.macresearch.org/no_need_to_parse_an_expression

The form for the shell script was based on examples at:
http://developer.apple.com/technotes/tn2002/tn2065.html
http://developer.apple.com/documentation/OpenSource/Conceptual/ 
ShellScripting/RegularExpressionsUnfettered/chapter_6_section_9.html
       	  l     ������  ��   	  
  
 i         I     ������
�� .aevtoappnull  �   � ****��  ��    k    ?       l     �� ��    > 8--------------------------------------------------------         l     �� ��    K EIn the "set _maxcontent" statement below, specify how many characters         l     �� ��    D >from the message to add to comments. Zero is a possible value.         l     �� ��    > 8--------------------------------------------------------         r         c          m     ���� <   m    ��
�� 
long  o      ���� 0 _maxcontent     ! " ! O     # $ # I  
 ������
�� .miscactvnull��� ��� null��  ��   $ m     % %Fnull     ߀��   1Yojimbo.app��@   �
  �����T �    �R�����e���R�@   @   ɥ   Cope   alis    �  Users                          NJ ����Yojimbo.app                                                    ����            ����  �cu             Applications krb  \/:Network:Servers:spitfire.genomecenter.ucdavis.edu:Users:keith:Applications krb:Yojimbo.app    Y o j i m b o . a p p    U s e r s  #/keith/Applications krb/Yojimbo.app   8/Network/Servers/spitfire.genomecenter.ucdavis.edu/Users 	 E Ecrbm  posx/Network/Servers/spitfire.genomecenter.ucdavis.edu/Users  ��   "  &�� & O   ? ' ( ' k   > ) )  * + * I   ������
�� .miscactvnull��� ��� null��  ��   +  , - , Z    ; . /���� . H    $ 0 0 l   # 1�� 1 I   #�� 2��
�� .coredoexbool        obj  2 4    �� 3
�� 
mvwr 3 m    ���� ��  ��   / k   ' 7 4 4  5 6 5 I  ' 4�� 7 8
�� .sysodlogaskr        TEXT 7 b   ' * 9 : 9 l 	 ' ( ;�� ; l 	 ' ( <�� < m   ' ( = = ; 5No Viewer Window is Open. Open one & select messages.   ��  ��   : o   ( )��
�� 
ret  8 �� > ?
�� 
appr > l 	 + , @�� @ m   + , A A  	Error....   ��   ? �� B C
�� 
givu B m   - .����  C �� D��
�� 
disp D m   / 0����  ��   6  E�� E L   5 7����  ��  ��  ��   -  F�� F O   <> G H G k   C= I I  J K J Z   C l L M���� L H   C M N N l  C L O�� O I  C L�� P��
�� .coredoexbool        obj  P 1   C H��
�� 
smgs��  ��   M k   P h Q Q  R S R I  P e�� T U
�� .sysodlogaskr        TEXT T b   P Y V W V b   P U X Y X l 	 P S Z�� Z m   P S [ [  No message was selected   ��   Y o   S T��
�� 
ret  W m   U X \ \  Transfer cancelled.    U �� ] ^
�� 
appr ] l 	 Z ] _�� _ m   Z ] ` `  Error...   ��   ^ �� a b
�� 
givu a m   ^ _����  b �� c��
�� 
disp c m   ` a����  ��   S  d�� d L   f h����  ��  ��  ��   K  e f e r   m v g h g 1   m r��
�� 
smgs h o      ���� 0 messagelist messageList f  i j i l  w w�� k��   k  Now loop through messages    j  l�� l X   w= m�� n m k   �8 o o  p q p l  � ��� r��   r  
Initialize    q  s t s r   � � u v u m   � � w w       v o      ���� 	0 _type   t  x y x r   � � z { z m   � � | |       { o      ���� 0 	_contents   y  } ~ } r   � �  �  m   � � � �       � o      ���� 	0 _name   ~  � � � r   � � � � � m   � � � �       � o      ���� 0 _person   �  � � � r   � � � � � m   � � � �       � o      ���� 0 _tag   �  � � � r   � � � � � m   � � � �       � o      ���� 	0 _body   �  � � � l  � ��� ���   � . (Classify message as incoming or outgoing    �  � � � r   � � � � � c   � � � � � n  � � � � � I   � ��� ����� 0 	countsent 	CountSent �  ��� � o   � ����� 0 m  ��  ��   �  f   � � � m   � ���
�� 
long � o      ���� 
0 _count   �  � � � Z   � � � ��� � � =  � � � � � o   � ����� 
0 _count   � m   � �����  � r   � � � � � m   � � � �  outgoing    � o      ���� 	0 _type  ��   � r   � � � � � m   � � � �  incoming    � o      ���� 	0 _type   �  � � � r   � � � � � n   � � � � � 1   � ���
�� 
subj � o   � ����� 0 m   � o      ���� 	0 _name   �  � � � r   � � � � � b   � � � � � o   � ����� 0 	_contents   � n  � � � � � I   � ��� ����� *0 generatemessagetext generateMessageText �  ��� � o   � ����� 0 m  ��  ��   �  f   � � � o      ���� 0 	_contents   �  � � � l  � ��� ���   �  Set up date sent    �  � � � r   � � � � n   � � � � 1   ���
�� 
drcv � o   � ����� 0 m   � o      ���� 	0 _date   �  � � � r   � � � c   � � � n   � � � 1  
��
�� 
year � o  
���� 	0 _date   � m  ��
�� 
TEXT � o      ���� 0 _yr   �  � � � r  & � � � c  " � � � n   � � � m  ��
�� 
mnth � o  ���� 	0 _date   � m  !��
�� 
nmbr � o      ���� 
0 _month   �  � � � l ''�� ���   �   Force month to format 'mm'    �  � � � r  '2 � � � c  '. � � � o  '*���� 
0 _month   � m  *-��
�� 
TEXT � o      ���� 	0 _smon   �  � � � Z  3L � ����� � A  3: � � � o  36���� 
0 _month   � m  69���� 
 � r  =H � � � b  =D � � � m  =@ � �  0    � o  @C���� 	0 _smon   � o      ���� 	0 _smon  ��  ��   �  � � � l MM�� ���   �  Force day to format 'dd'    �  � � � r  M\ � � � c  MX � � � n  MT � � � 1  PT��
�� 
day  � o  MP���� 	0 _date   � m  TW��
�� 
TEXT � o      ���� 0 _day   �  � � � Z  ]z � ����� � A  ]h � � � n  ]d � � � 1  `d��
�� 
day  � o  ]`���� 	0 _date   � m  dg���� 
 � r  kv � � � b  kr � � � m  kn � �  0    � o  nq���� 0 _day   � o      ���� 0 _day  ��  ��   �    l {{����   . ( Create date sent with yyyy/mm/dd format     r  {� b  {� b  {�	
	 b  {� b  {� o  {~�� 0 _yr   m  ~�  /    o  ���~�~ 	0 _smon  
 m  ��  /    o  ���}�} 0 _day   o      �|�| 
0 _datef    l ���{�{   " Get sender's name or address     r  �� I ���z�y
�z .emaleafnTEXT        TEXT n  �� 1  ���x
�x 
sndr o  ���w�w 0 m  �y   o      �v�v 0 _sender    l ���u�u   : 4If there is no name, _sender defaults to the address     l ���t �t    = 7If the message was received, set the name to the sender    !"! l ���s#�s  # %  and set a tag variable to "in"   " $%$ Z  ��&'�r�q& = ��()( o  ���p�p 	0 _type  ) m  ��**  incoming   ' k  ��++ ,-, r  ��./. b  ��010 m  ��22 	 <:    1 o  ���o�o 0 _sender  / o      �n�n 0 _person  - 3�m3 r  ��454 c  ��676 m  ��88  in   7 m  ���l
�l 
TEXT5 o      �k�k 0 _tag  �m  �r  �q  % 9:9 l ���j;�j  ; > 8If the message was sent, get the first recipient's name,   : <=< l ���i>�i  > ) #if available, otherwise the address   = ?@? Z  �)AB�h�gA = ��CDC o  ���f�f 	0 _type  D m  ��EE  outgoing   B k  �%FF GHG r  ��IJI n  ��KLK 2 ���e
�e 
trcpL o  ���d�d 0 m  J o      �c�c 0 therecipients theRecipientsH MNM Z  �OP�bQO l ��R�aR I ���`S�_
�` .coredoexbool        obj S n  ��TUT 1  ���^
�^ 
pnamU n  ��VWV 4  ���]X
�] 
cobjX m  ���\�\ W o  ���[�[ 0 therecipients theRecipients�_  �a  P r  �YZY b  ��[\[ m  ��]] 	 >:    \ n  ��^_^ 1  ���Z
�Z 
pnam_ n  ��`a` 4  ���Yb
�Y 
cobjb m  ���X�X a o  ���W�W 0 therecipients theRecipientsZ o      �V�V 0 _person  �b  Q r  cdc b  efe m  gg 	 >:    f n  hih 1  �U
�U 
raddi n  jkj 4  �Tl
�T 
cobjl m  �S�S k o  �R�R 0 therecipients theRecipientsd o      �Q�Q 0 _person  N m�Pm r  %non c  !pqp m  rr 	 out   q m   �O
�O 
TEXTo o      �N�N 0 _tag  �P  �h  �g  @ sts l **�Mu�M  u ) #Get message body to add to contents   t vwv r  *1xyx m  *-zz      y o      �L�L 	0 _body  w {|{ r  2C}~} c  2?� b  2;��� o  25�K�K 	0 _body  � n  5:��� 1  6:�J
�J 
ctnt� o  56�I�I 0 m  � m  ;>�H
�H 
ctxt~ o      �G�G 	0 _body  | ��� r  DO��� n  DK��� 1  GK�F
�F 
leng� o  DG�E�E 	0 _body  � o      �D�D 0 _length  � ��� r  PW��� m  PS��      � o      �C�C 0 
_shorttext  � ��� Z  X����B�A� ?  X[��� o  XY�@�@ 0 _maxcontent  � m  YZ�?�?  � k  ^��� ��� r  ^o��� n ^k��� I  _k�>��=�> 0 min  � ��� o  _b�<�< 0 _length  � ��;� [  bg��� o  bc�:�: 0 _maxcontent  � m  cf�9�9 �;  �=  �  f  ^_� o      �8�8 
0 _allow  � ��� l pp�7��7  � , &allow for spaces that will be replaced   � ��� r  p���� e  p��� n  p���� 7 s��6��
�6 
ctxt� m  y{�5�5 � o  |��4�4 
0 _allow  � o  ps�3�3 	0 _body  � o      �2�2 0 _toomuch  � ��� l ���1�0�1  �0  � ��� l ���/��/  � = 7 Set shell script to remove tabs, returns, extra spaces   � ��� l ���.��.  � A ; **********************************************************   � ��� r  ����� I ���-��,
�- .sysoexecTEXT���     TEXT� b  ����� b  ����� b  ����� b  ����� b  ����� l 	����+� m  ����  echo    �+  � n  ����� 1  ���*
�* 
strq� o  ���)�) 0 _toomuch  � l 	����(� l 	����'� m  ���� D >| perl -e "while(\$line = <STDIN>) {-e \$line =~ s/[	
]+/ /g ;   �'  �(  � m  ����  print \$line; }"   � l 	����&� l 	����%� m  ���� B <| perl -e "while(\$line = <STDIN>) {-e \$line =~ s/ +/ /g ;    �%  �&  � m  ����  print \$line; }"   �,  � o      �$�$ 
0 _text1  � ��� l ���#��#  � A ;-**********************************************************   � ��� l ���"�!�"  �!  � ��� l ��� ��   � < 6Now find the maximum length of the content to be added   � ��� r  ����� n  ����� 1  ���
� 
leng� o  ���� 
0 _text1  � o      �� 0 
_lengthnow  � ��� r  ����� n ����� I  ������ 0 min  � ��� o  ���� 0 
_lengthnow  � ��� o  ���� 0 _maxcontent  �  �  �  f  ��� o      �� 	0 _most  � ��� r  ����� e  ���� n  ����� 7�����
� 
ctxt� m  ���� � o  ���� 	0 _most  � o  ���� 
0 _text1  � o      �� 0 
_shorttext  � ��� r  ����� b  ����� m  ���� 	 ---   � o  ���� 0 
_shorttext  � o      �� 0 
_shorttext  �  �B  �A  � ��� l �����  � 9 3-Now set up the note with body, comments, and tags.   � ��� O  �8��� k  �7��    r  �$ I � ��
� .corecrel****      � null�   �

�
 
kocl m  ���	
�	 
YNot ��
� 
prdt l 	�� l 
�	�	 K  �

 �
� 
pcnt o  � �� 0 	_contents   �
� 
pnam o  �� 	0 _name   � ��
�  
Cmnt b  	 b  	 b  	 o  	���� 
0 _datef   m         o  ���� 0 _person   o  ���� 0 
_shorttext  ��  �  �  �   o      ���� 0 anewnote aNewNote �� I %7��
�� .CopeAddTnull���    cct6 J  %-  m  %(  email    �� o  (+���� 0 _tag  ��   ����
�� 
to   o  03���� 0 anewnote aNewNote��  ��  � m  �� %�  �� 0 m   n o   z }���� 0 messagelist messageList��   H 4   < @�� 
�� 
mvwr  m   > ?���� ��   ( m    !!�null     ߀��  <Mail.app.app�@   �
  �����T �    �R�����e���R�@   @   ɥ   emal   alis    2  raiden                     �%��H+    <Mail.app                                                         ��L��        ����  	                Applications    �&6      �M*`      <  raiden:Applications:Mail.app    M a i l . a p p    r a i d e n  Applications/Mail.app   / ��  ��    "#" l     ������  ��  # $%$ l     ��&��  &   Function routines follow   % '(' l     ������  ��  ( )*) l     ��+��  + ; 5Set up a function to generate the body of the message   * ,-, l     ��.��  .  (by Jim Correia)   - /0/ i    121 I      ��3���� *0 generatemessagetext generateMessageText3 4��4 o      ���� 0 m  ��  ��  2 O     D565 k    C77 898 r    	:;: n    <=< 1    ��
�� 
sndr= o    ���� 0 m  ; o      ���� 0 _sender  9 >?> r   
 @A@ n   
 BCB 1    ��
�� 
subjC o   
 ���� 0 m  A o      ���� 0 _subject  ? DED r    FGF c    HIH n    JKJ 1    ��
�� 
rdrcK o    ���� 0 m  I m    ��
�� 
TEXTG o      ���� 	0 _date  E LML r    NON n    PQP 1    ��
�� 
ctntQ o    ���� 0 m  O o      ���� 0 	_contents  M RSR r    %TUT b    #VWV b    !XYX m    ZZ  From:    Y o     ���� 0 _sender  W o   ! "��
�� 
ret U o      ����  0 _messagestring _messageStringS [\[ r   & /]^] b   & -_`_ b   & +aba b   & )cdc o   & '����  0 _messagestring _messageStringd m   ' (ee  	Subject:    b o   ) *���� 0 _subject  ` o   + ,��
�� 
ret ^ o      ����  0 _messagestring _messageString\ fgf r   0 9hih b   0 7jkj b   0 5lml b   0 3non o   0 1����  0 _messagestring _messageStringo m   1 2pp  Date:    m o   3 4���� 	0 _date  k o   5 6��
�� 
ret i o      ����  0 _messagestring _messageStringg q��q r   : Crsr b   : Atut b   : ?vwv b   : =xyx o   : ;����  0 _messagestring _messageStringy o   ; <��
�� 
ret w o   = >��
�� 
ret u o   ? @���� 0 	_contents  s o      ����  0 _messagestring _messageString��  6 m     !0 z{z l     ������  ��  { |}| l     ��~��  ~ 1 +Function to pass a variable which indicates   } � l     �����  � ' !whether I sent the message or not   � ��� l     ������  ��  � ��� i    ��� I      ������� 0 	countsent 	CountSent� ���� o      ���� 0 m  ��  ��  � k     D�� ��� O     A��� k    @�� ��� r    ��� I   �����
�� .emaleauaTEXT        TEXT� n    ��� 1    ��
�� 
sndr� o    ���� 0 m  ��  � o      ���� 
0 _saddr  � ��� r    ��� c    ��� m    ����  � m    ��
�� 
long� o      ���� 0 _ctsent  � ���� X    @����� Z   & ;������� E   & +��� n   & )��� 1   ' )��
�� 
emad� o   & '���� 	0 _acct  � o   ) *���� 
0 _saddr  � k   . 7�� ��� r   . 5��� c   . 3��� [   . 1��� o   . /���� 0 _ctsent  � m   / 0���� � m   1 2��
�� 
long� o      ���� 0 _ctsent  � ��� l  6 6�����  � 6 0 _ctsent counts the number of times the sender's   � ���� l  6 6�����  � + % address matches one of mine (0 or 1)   ��  ��  ��  �� 	0 _acct  � 2   ��
�� 
mact��  � m     !� ���� L   B D�� o   B C���� 0 _ctsent  ��  � ��� l     ������  ��  � ��� i    ��� I      ������� 0 min  � ��� o      ���� 0 x  � ���� o      ���� 0 y  ��  ��  � Z     ������ l    ���� B    ��� o     ���� 0 x  � o    ���� 0 y  ��  � L    �� o    ���� 0 x  ��  � L    �� o    ���� 0 y  � ���� l     ������  ��  ��       ���������  � ��������
�� .aevtoappnull  �   � ****�� *0 generatemessagetext generateMessageText�� 0 	countsent 	CountSent�� 0 min  � �� ��������
�� .aevtoappnull  �   � ****��  ��  � ���� 0 m  � j������ %��!���� =���� A��������~ [ \ `�}�|�{�z w�y |�x ��w ��v ��u ��t�s�r � ��q�p�o�n�m�l�k�j�i�h�g�f ��e�d ��c�b�a�`*28E�_�^�]]g�\rz�[�Z�Y�X��W�V�U�T�S��R�����Q�P�O�N��M�L�K�J�I�H�G�F�E�� <
�� 
long�� 0 _maxcontent  
�� .miscactvnull��� ��� null
�� 
mvwr
�� .coredoexbool        obj 
�� 
ret 
�� 
appr
�� 
givu
�� 
disp�� 
� .sysodlogaskr        TEXT
�~ 
smgs�} 0 messagelist messageList
�| 
kocl
�{ 
cobj
�z .corecnte****       ****�y 	0 _type  �x 0 	_contents  �w 	0 _name  �v 0 _person  �u 0 _tag  �t 	0 _body  �s 0 	countsent 	CountSent�r 
0 _count  
�q 
subj�p *0 generatemessagetext generateMessageText
�o 
drcv�n 	0 _date  
�m 
year
�l 
TEXT�k 0 _yr  
�j 
mnth
�i 
nmbr�h 
0 _month  �g 	0 _smon  �f 

�e 
day �d 0 _day  �c 
0 _datef  
�b 
sndr
�a .emaleafnTEXT        TEXT�` 0 _sender  
�_ 
trcp�^ 0 therecipients theRecipients
�] 
pnam
�\ 
radd
�[ 
ctnt
�Z 
ctxt
�Y 
leng�X 0 _length  �W 0 
_shorttext  �V �U 0 min  �T 
0 _allow  �S 0 _toomuch  
�R 
strq
�Q .sysoexecTEXT���     TEXT�P 
0 _text1  �O 0 
_lengthnow  �N 	0 _most  
�M 
YNot
�L 
prdt
�K 
pcnt
�J 
Cmnt�I 
�H .corecrel****      � null�G 0 anewnote aNewNote
�F 
to  
�E .CopeAddTnull���    cct6��@��&E�O� *j UO�+*j O*�k/j  ��%���m�j� OhY hO*�k/�*a ,j  a �%a %�a �m�j� OhY hO*a ,E` O�_ [a a l kh  a E` Oa E` Oa E` Oa E` Oa  E` !Oa "E` #O)�k+ $�&E` %O_ %k  a &E` Y 	a 'E` O�a (,E` O_ )�k+ )%E` O�a *,E` +O_ +a ,,a -&E` .O_ +a /,a 0&E` 1O_ 1a -&E` 2O_ 1a 3 a 4_ 2%E` 2Y hO_ +a 5,a -&E` 6O_ +a 5,a 3 a 7_ 6%E` 6Y hO_ .a 8%_ 2%a 9%_ 6%E` :O�a ;,j <E` =O_ a >  a ?_ =%E` Oa @a -&E` !Y hO_ a A  Y�a B-E` CO_ Ca k/a D,j  a E_ Ca k/a D,%E` Y a F_ Ca k/a G,%E` Oa Ha -&E` !Y hOa IE` #O_ #�a J,%a K&E` #O_ #a L,E` MOa NE` OO�j �)_ M�a Pl+ QE` RO_ #[a K\[Zk\Z_ R2EE` SOa T_ Sa U,%a V%a W%a X%a Y%j ZE` [O_ [a L,E` \O)_ \�l+ QE` ]O_ [[a K\[Zk\Z_ ]2EE` OOa ^_ O%E` OY hO� I*a a _a `a a_ a D_ a b_ :a c%_ %_ O%�a d eE` fOa g_ !lva h_ fl iU[OY�OUU� �D2�C�B���A�D *0 generatemessagetext generateMessageText�C �@��@ �  �?�? 0 m  �B  � �>�=�<�;�:�9�> 0 m  �= 0 _sender  �< 0 _subject  �; 	0 _date  �: 0 	_contents  �9  0 _messagestring _messageString� 
!�8�7�6�5�4Z�3ep
�8 
sndr
�7 
subj
�6 
rdrc
�5 
TEXT
�4 
ctnt
�3 
ret �A E� A��,E�O��,E�O��,�&E�O��,E�O�%�%E�O��%�%�%E�O��%�%�%E�O��%�%�%E�U� �2��1�0���/�2 0 	countsent 	CountSent�1 �.��. �  �-�- 0 m  �0  � �,�+�*�)�, 0 m  �+ 
0 _saddr  �* 0 _ctsent  �) 	0 _acct  � 	!�(�'�&�%�$�#�"�!
�( 
sndr
�' .emaleauaTEXT        TEXT
�& 
long
�% 
mact
�$ 
kocl
�# 
cobj
�" .corecnte****       ****
�! 
emad�/ E� >��,j E�Oj�&E�O +*�-[��l kh ��,� �k�&E�OPY h[OY��UO�� � �������  0 min  � ��� �  ��� 0 x  � 0 y  �  � ��� 0 x  � 0 y  �  � �� �Y � ascr  ��ޭ